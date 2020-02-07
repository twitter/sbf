/**
 * Copyright 2020 Twitter. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.twitter.sbf.core;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.BiFunction;
import java.util.function.Function;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.math3.distribution.GumbelDistribution;
import org.apache.commons.math3.random.RandomAdaptor;

import com.twitter.sbf.graph.ConnectedComponents;
import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.graph.SubGraphWithoutWeakLinks;
import com.twitter.sbf.util.Range;
import com.twitter.sbf.util.Util;

import it.unimi.dsi.fastutil.PriorityQueue;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;
import it.unimi.dsi.fastutil.objects.ObjectList;


/**
 * Class that computes sparse binary matrix factorization of an input graph using
 * Metropolis-Hastings sampling algorithm (can also be viewed as Simulated Annealing).
 */
public class MHAlgorithm {

  private Graph g;
  private AlgorithmConfig config;
  private SparseBinaryMatrix initMatrix;
  private PrintWriter diagnosticsWriter;
  private long emptyProposals = 0;
  private long sameProposals = 0;
  private GumbelDistribution gumbel;
  private Lock[] colLocks;

  public MHAlgorithm(AlgorithmConfig config,
                     Graph g,
                     SparseBinaryMatrix matrix,
                     PrintWriter diagnosticsWriter) {
    this.g = g;
    this.config = config;
    this.initMatrix = matrix;
    this.diagnosticsWriter = diagnosticsWriter;
    gumbel = new GumbelDistribution(config.rng, 0, 1);
    if (config.updateImmediately && !config.noLocking) {
      colLocks = new ReentrantLock[config.k];
      for (int i = 0; i < colLocks.length; i++) {
        colLocks[i] = new ReentrantLock();
      }
    }
  }

  /**
   * Optimize initial matrix using Metropolis-Hastings sampling algorithm
   *
   * @return optimized matrix
   */
  public SparseBinaryMatrix optimize() {
    double wtCoeff = config.wtCoeff;
    diagnosticsWriter.println("WeightCoeff: " + wtCoeff);

    diagnosticsWriter.println(Util.repeatString(100, "-"));
    diagnosticsWriter.println(headerLine(metricFieldsToUse(g, config)));
    diagnosticsWriter.println(
        getMetricsLine(initMatrix, 0, 0, 0, 0, true)
    );

    int epoch;
    for (epoch = 1; epoch <= config.maxEpoch; epoch++) {
      int acceptCount = runEpoch(initMatrix, epoch,
        config.maxEpoch < 1000 || epoch % config.evalEvery == 0);
      double acceptRate = acceptCount * 1.0 / g.getNumVertices();
      if (!config.runAllEpochs && acceptRate < config.eps) {
        break;
      }
    }

    if (config.divideResultIntoConnectedComponents) {
      SparseBinaryMatrix dividedZ = divideIntoConnectedComponents(initMatrix);
      // divideIntoConnectedComponents already removes clusters smaller than 2 i.e. singletons
      if (config.minClusterSize > 2) {
        return removeSmallClusters(dividedZ, config.minClusterSize);
      } else {
        return dividedZ;
      }
    } else {
      if (config.minClusterSize > 1) {
        return removeSmallClusters(initMatrix, config.minClusterSize);
      } else {
        return initMatrix;
      }
    }
  }

  private SparseBinaryMatrix removeSmallClusters(SparseBinaryMatrix z, int minClusterSize) {
    diagnosticsWriter.format(
        "Going to remove columns with fewer than %d nonzeros\n", minClusterSize
    );
    diagnosticsWriter.flush();

    int numRemoved = 0;
    ImmutableList.Builder<IntSet> newCols = new ImmutableList.Builder<>();
    for (int i = 0; i < z.getNumCols(); i++) {
      if (z.getColumn(i).size() >= minClusterSize) {
        newCols.add(z.getColumn(i));
      } else {
        numRemoved++;
      }
    }

    ImmutableList<IntSet> newColsList = newCols.build();
    IntSet[] newColsArray =
        Arrays.copyOf(newColsList.toArray(), newColsList.size(), IntSet[].class);
    SparseBinaryMatrix newZ = new SparseBinaryMatrix(g.getNumVertices(), newColsArray.length);
    newZ.initFromColSets(newColsArray);

    diagnosticsWriter.println(getMetricsLine(newZ, 0, 0, 0, 0, true));
    diagnosticsWriter.format("Removed %d clusters which are smaller than %d\n",
        numRemoved, minClusterSize);
    diagnosticsWriter.flush();

    return newZ;
  }

  private SparseBinaryMatrix divideIntoConnectedComponents(SparseBinaryMatrix z) {
    diagnosticsWriter.print("Going to divide matrix into connected components");
    if (config.removeWeakLinksForConnectedComponents) {
      diagnosticsWriter.println(" (removing weak links)");
    } else {
      diagnosticsWriter.println();
    }
    diagnosticsWriter.flush();

    ImmutableList.Builder<IntSet> newCols = new ImmutableList.Builder<>();
    long tic = System.nanoTime();
    int numOriginalSingletons = 0;
    int numNewSingletons = 0;
    int numEmpty = 0;
    long numWeakLinks = 0;
    for (int i = 0; i < z.getNumCols(); i++) {
      if (z.getColumn(i).size() > 1) {
        Set<IntSet> comps;
        if (config.removeWeakLinksForConnectedComponents) {
          SubGraphWithoutWeakLinks sg = new SubGraphWithoutWeakLinks(g, z.getColumn(i));
          comps = ConnectedComponents.connectedComponents(sg);
          numWeakLinks += sg.getNumWeakLinks();
        } else {
          comps = ConnectedComponents.subGraphConnectedComponents(g, z.getColumn(i));
        }

        for (IntSet c : comps) {
          newCols.add(c);
          if (c.size() == 1) {
            numNewSingletons += 1;
          }
        }
      } else if (z.getColumn(i).size() == 1) {
        numOriginalSingletons++;
      } else {
        numEmpty++;
      }
    }
    ImmutableList<IntSet> newColsList = newCols.build();
    IntSet[] newColsArray =
        Arrays.copyOf(newColsList.toArray(), newColsList.size(), IntSet[].class);
    SparseBinaryMatrix dividedZ = new SparseBinaryMatrix(g.getNumVertices(), newColsArray.length);
    dividedZ.initFromColSets(newColsArray);
    long toc = System.nanoTime();
    double updateTime = (toc - tic) * 1e-9;
    if (config.removeWeakLinksForConnectedComponents) {
      diagnosticsWriter.println("Number of weak links found: " + numWeakLinks);
    }
    diagnosticsWriter.println(
        getMetricsLine(dividedZ, 0, updateTime, 0, 0,
            true)
    );

    int originalEmptyRows = (int) Math.round(z.emptyRowProportion() * z.getNumRows());
    int newEmptyRows = (int) Math.round(dividedZ.emptyRowProportion() * dividedZ.getNumRows());
    diagnosticsWriter.format(
        "Found %d original singletons, %d new singletons, %d empty columns\n",
        numOriginalSingletons, numNewSingletons, numEmpty);
    diagnosticsWriter.format("Dividing into connected components and removing singletons "
        + "and empties changes the number of empty rows (unassigned nodes) from %d to %d\n",
        originalEmptyRows, newEmptyRows);
    diagnosticsWriter.format(
        "Dividing into connected components and removing singletons "
            + "and empties changes K from %d to %d\n",
        z.getNumCols(), dividedZ.getNumCols());
    diagnosticsWriter.flush();
    return dividedZ;
  }

  private Integer[] shuffleVertexIds() {
    Integer[] shuffledIds = new Integer[g.getNumVertices()];
    for (int i = 0; i < g.getNumVertices(); i++) {
      shuffledIds[i] = i;
    }
    Collections.shuffle(Arrays.asList(shuffledIds));
    return shuffledIds;
  }

  private IntSet lockAndReturnLockedColumns(SparseBinaryMatrix matrix, int vertexId) {
    IntSet colsToLock = new IntOpenHashSet();
    synchronized (this) {
      for (int nId : g.getNeighbors(vertexId)) {
        // SUPPRESS CHECKSTYLE NestedForDepth
        for (int colId : matrix.getRow(nId)) {
          if (!colsToLock.contains(colId)) {
            colsToLock.add(colId);
            colLocks[colId].lock();
          }
        }
      }
    }
    return colsToLock;
  }

  private int runEpochSerial(
      SparseBinaryMatrix matrix,
      Integer[] vertexIds,
      double multiplier,
      double wtCoeff) {
    int acceptCount = 0;
    for (int vertexId : vertexIds) {
      Optional<int[]> newRow = mhStep(matrix, vertexId, multiplier, wtCoeff);
      if (newRow.isPresent()) {
        acceptCount++;
        matrix.updateRow(vertexId, newRow.get());
      }
    }
    return acceptCount;
  }

  private int runEpoch(SparseBinaryMatrix matrix, final int epoch, boolean doMetrics) {
    double wtCoeff = config.wtCoeff;

    long tic = System.nanoTime();
    emptyProposals = 0;
    sameProposals = 0;
    double temp = Math.pow(config.temperatureRatio, epoch) * config.maxTemperature;
    double multiplier = config.useTemperatureSchedule ? 1.0 / temp : config.scaleCoeff;

    int n = g.getNumVertices();

    Integer[] shuffledIds = shuffleVertexIds();

    final Hashtable<Integer, Optional<int[]>> updates = new Hashtable<>(n, 1.0f);

    int acceptCount = 0;
    if (config.cpu <= 1) {
      acceptCount = runEpochSerial(matrix, shuffledIds, multiplier, wtCoeff);
    } else {
      /*
       * All of the below complicated logic is complicated because it's multi-threaded and there's
       * a few different options to enable/disable when running multi-threaded;
       * the essential logic is the same as in runEpochSerial.
       */
      long maxWork = 0;
      long minWork = Long.MAX_VALUE;
      int numTasks = Math.max(1, n / 1000);
      ArrayList<Future<Integer>> futures = new ArrayList<>(numTasks);
      ExecutorService exec = Executors.newWorkStealingPool(config.cpu);
      for (Range range : Util.chunkArray(n, numTasks)) {
        long workHere = 0;
        for (int i = range.start; i < range.end; i++) {
          workHere += g.getDegree(shuffledIds[i]);
        }
        if (workHere > maxWork) {
          maxWork = workHere;
        }
        if (workHere < minWork && workHere > 0) {
          minWork = workHere;
        }
        Future<Integer> fut = exec.submit(() -> {
          int accepted = 0;
          for (int i = range.start; i < range.end; i++) {
            int vertexId = shuffledIds[i];
            if (!config.updateImmediately) {
              Optional<int[]> newRow = mhStep(matrix, vertexId, multiplier, wtCoeff);
              if (newRow.isPresent()) {
                accepted++;
              }
              updates.put(vertexId, newRow);
            } else {
              if (!config.noLocking) {
                IntSet colsToLock = lockAndReturnLockedColumns(matrix, vertexId);
                try {
                  Optional<int[]> newRow =
                      mhStep(matrix, vertexId, multiplier, wtCoeff);
                  if (newRow.isPresent()) {
                    matrix.updateRow(vertexId, newRow.get());
                    accepted++;
                  }
                } finally {
                  for (int colId : colsToLock) {
                    colLocks[colId].unlock();
                  }
                }
              } else { // HogWild!!
                Optional<int[]> newRow = mhStep(matrix, vertexId, multiplier, wtCoeff);
                if (newRow.isPresent()) {
                  matrix.updateRow(vertexId, newRow.get());
                  accepted++;
                }
              }
            }
          }
          return accepted;
        });
        futures.add(fut);
      }

      if (maxWork / minWork >= 2) {
        diagnosticsWriter.format("epoch is %d, maxWork is %d, minWork is %d\n",
            epoch, maxWork, minWork);
        diagnosticsWriter.flush();
      }

      exec.shutdown();
      try {
        while (!exec.isTerminated()) {
          exec.awaitTermination(1200, TimeUnit.SECONDS);
          if (!exec.isTerminated()) {
            diagnosticsWriter.println("epoch is " + epoch + ", going to wait for 20 more minutes, "
                + "executor.isTerminated " + exec.isTerminated());
            diagnosticsWriter.flush();
          }
        }
        for (Future<Integer> fut : futures) {
          acceptCount += fut.get();
        }
      } catch (InterruptedException | ExecutionException e) {
        e.printStackTrace(diagnosticsWriter);
      }

      if (!config.updateImmediately) {
        for (Map.Entry<Integer, Optional<int[]>> update : updates.entrySet()) {
          if (update.getValue().isPresent()) {
            matrix.updateRow(update.getKey(), update.getValue().get());
          }
        }
      }
    }

    if (doMetrics) {
      double acceptRate = acceptCount * 1.0 / n;

      long toc = System.nanoTime();
      double timeInSecs = (toc - tic) * 1e-9;
      diagnosticsWriter.println(
          getMetricsLine(matrix, epoch, timeInSecs, acceptRate, 1.0 / multiplier,
              epoch % config.evalEvery == 0)
      );
      diagnosticsWriter.flush();
    }

    return acceptCount;
  }

  private static final Map<String, Integer> METRICS_FIELDS_TO_WIDTHS =
      ImmutableMap.<String, Integer>builder()
          .put("epoch", 6)
          .put("nnz/vertex", 12)
          .put("emptyRow", 9)
          .put("emptyCol", 9)
          .put("orphans", 7)
          .put("minSize", 8)
          .put("maxSize", 8)
          .put("prec", 5)
          .put("rec", 5)
          .put("f1", 5)
          .put("wtPrec", 6)
          .put("wtRec", 6)
          .put("wtF1", 6)
          .put("error", 6)
          .put("totalSec", 9)
          .put("updateSec", 9)
          .put("evalSec", 8)
          .put("acceptRate", 10)
          .put("temp", 6)
          .put("emptyProposal", 13)
          .put("sameProposal", 13)
          .build();

  private static final Map<String, String> METRICS_FIELDS_TO_FORMATS =
      ImmutableMap.<String, String>builder()
          .put("epoch", "%6d")
          .put("nnz/vertex", "%12.2g")
          .put("emptyRow", "%9.2g")
          .put("emptyCol", "%9.2g")
          .put("orphans", "%7.2f")
          .put("minSize", "%8d")
          .put("maxSize", "%8d")
          .put("prec", "%5.1f")
          .put("rec", "%5.1f")
          .put("f1", "%5.1f")
          .put("wtPrec", "%6.1f")
          .put("wtRec", "%6.1f")
          .put("wtF1", "%6.1f")
          .put("error", "%6.2f")
          .put("totalSec", "%9.2f")
          .put("updateSec", "%9.2f")
          .put("evalSec", "%8.2f")
          .put("acceptRate", "%10.4f")
          .put("temp", "%6.4g")
          .put("emptyProposal", "%13.2f")
          .put("sameProposal", "%13.2f")
          .build();

  private static final List<String> METRICS_FIELDS = Lists.newArrayList(
      "epoch",
      "nnz/vertex",
      "emptyRow",
      "orphans",
      "emptyCol",
      "minSize",
      "maxSize",
      "prec",
      "rec",
      "f1",
      "wtPrec",
      "wtRec",
      "wtF1",
      //"error",
      "totalSec",
      //"updateSec",
      "evalSec",
      "acceptRate",
      "temp",
      "emptyProposal",
      "sameProposal"
  );

  private Map<String, Object> getMetrics(SparseBinaryMatrix matrix,
                                         int epoch,
                                         double elapsedSeconds,
                                         double acceptRate,
                                         double temperature,
                                         boolean doExpensiveEval) {
    Graph graph = g;
    ImmutableMap.Builder<String, Object> metricsMap = new ImmutableMap.Builder<>();
    metricsMap.put("epoch", epoch);
    metricsMap.put("nnz/vertex", matrix.nnzPerRow());
    metricsMap.put("emptyRow", matrix.emptyRowProportion());
    metricsMap.put("emptyCol", matrix.emptyColProportion());
    metricsMap.put("minSize", matrix.minColSizeAboveZero());
    metricsMap.put("maxSize", matrix.maxColSize());
    metricsMap.put("totalSec", elapsedSeconds);
    metricsMap.put("acceptRate", acceptRate);
    metricsMap.put("temp", temperature);
    metricsMap.put("emptyProposal", emptyProposals * 1.0 / graph.getNumVertices());
    metricsMap.put("sameProposal", sameProposals * 1.0 / graph.getNumVertices());
    if (doExpensiveEval) {
      long tic = System.nanoTime();
      // Evaluate prediction metrics
      PredictionStat stat;
      if (graph.getNumEdges() < 5000L) {
        int subsetSize =
            (config.evalRatio < 1.0) ? ((int) (graph.getNumVertices() * config.evalRatio))
                : graph.getNumVertices();
        int[] evalVertexIds =
            (config.evalRatio < 1.0) ? Util.reservoirSampling(graph.getNumVertices(), subsetSize)
                : graph.getAllVertexIds();
        // Parallel execution using Executor
        ExecutorService executor = Executors.newFixedThreadPool(config.cpu);
        stat = getPredictionStatVertexSampling(graph, matrix, executor, evalVertexIds);
        executor.shutdown();
        try {
          while (!executor.isTerminated()) {
            executor.awaitTermination(20, TimeUnit.MINUTES);
          }
        } catch (InterruptedException e) {
          throw new RuntimeException(e);
        }
        metricsMap.put("orphans", stat.orphanRate());
        metricsMap.put("prec", stat.precision() * 100);
        metricsMap.put("rec", stat.recall() * 100);
        metricsMap.put("f1", stat.f1() * 100);
        metricsMap.put("wtPrec", stat.weightedPrecision() * 100);
        metricsMap.put("wtRec", stat.weightedRecall() * 100);
        metricsMap.put("wtF1", stat.weightedF1() * 100);
      } else {
        Map<String, PredictionStat> ss = getPredictionStatEdgeSampling(graph, matrix, config.rng,
            100000, 1000000, 5000);
        double precision = ss.get("precision").precision();
        double recall = ss.get("recall").recall();
        double wtPrec = ss.get("precision").weightedPrecision();
        double wtRec = ss.get("recall").weightedRecall();
        metricsMap.put("orphans", ss.get("orphans").orphanRate());
        metricsMap.put("prec", precision * 100);
        metricsMap.put("rec", recall * 100);
        metricsMap.put("f1", 200 * precision * recall / (precision + recall));
        metricsMap.put("wtPrec", wtPrec * 100);
        metricsMap.put("wtRec", wtRec * 100);
        metricsMap.put("wtF1", 200 * wtPrec * wtRec / (wtPrec + wtRec));
      }

      long toc = System.nanoTime();
      double evalTime = (toc - tic) * 1e-9;
      metricsMap.put("evalSec", evalTime);
    }
    return metricsMap.build();
  }

  private String getMetricsLine(SparseBinaryMatrix matrix,
                                int epoch,
                                double elapsedSeconds,
                                double acceptRate,
                                double temperature,
                                boolean doExpensiveEval) {
    return metricsLineFromValues(
        metricFieldsToUse(g, config),
        getMetrics(matrix, epoch, elapsedSeconds, acceptRate, temperature, doExpensiveEval)
    );
  }


  private static List<String> metricFieldsToUse(Graph graph, AlgorithmConfig config) {
    List<String> ret = METRICS_FIELDS;
    if (!config.useTemperatureSchedule) {
      ret.remove("temp");
    }
    if (!graph.isWeighted()) {
      ret.remove("wtPrec");
      ret.remove("wtRec");
      ret.remove("wtF1");
    }
    if (graph.getNumVertices() < 20000) {
      ret.remove("evalSec");
      ret.remove("updateSec");
    }
    return ret;
  }

  private static String metricsLineFromValues(
      List<String> fieldsToUse,
      Map<String, Object> valueMap) {
    ImmutableMap.Builder<String, String> builder = ImmutableMap.<String, String>builder();
    for (Map.Entry<String, Object> entry : valueMap.entrySet()) {
      builder.put(
          entry.getKey(),
          String.format(METRICS_FIELDS_TO_FORMATS.get(entry.getKey()), entry.getValue())
      );
    }
    return metricsLineFromStrings(fieldsToUse, builder.build());
  }

  private static String metricsLineFromStrings(
      List<String> fieldsToUse,
      Map<String, String> valueMap) {
    StringBuilder ret = new StringBuilder();
    for (String field : fieldsToUse) {
      int width = METRICS_FIELDS_TO_WIDTHS.get(field);
      String toPrint = valueMap.get(field);
      if (toPrint == null) {
        toPrint = "-";
      }
      ret.append(String.format("%" + width + "s ", toPrint));
    }
    return ret.toString();
  }

  private static String headerLine(List<String> fieldsToUse) {
    return metricsLineFromStrings(fieldsToUse, Maps.toMap(fieldsToUse, f -> f));
  }

  private static boolean isOrphan(Graph graph,
                                  SparseBinaryMatrix matrix,
                                  int vId,
                                  int[] vFactor) {
    boolean isOrphan = true;
    if (vFactor.length == 0) {
      isOrphan = false;
    }
    for (int neighbor : graph.getNeighbors(vId)) {
      if (Util.hasCommonElement(vFactor, matrix.getRow(neighbor))) {
        isOrphan = false;
        break;
      }
    }
    return isOrphan;
  }

  private static Map<String, PredictionStat> getPredictionStatEdgeSampling(
      Graph graph,
      SparseBinaryMatrix matrix,
      RandomAdaptor rng,
      int recallSamples,
      int precisionSamples,
      int verticesToSample) {
    PredictionStat recallStat = new PredictionStat();
    for (int i = 0; i < precisionSamples; i++) {
      // sample first vertex according to degree distribution
      int v1 = graph.getDegreeDistribution(rng).sample();
      if (graph.getDegree(v1) > 0) {
        int v2 = graph.getNeighbors(v1)[rng.nextInt(graph.getDegree(v1))];
        float edgeWeight = graph.getWeightOfEdge(v1, v2);
        boolean predictEdge = Util.hasCommonElement(matrix.getRow(v1), matrix.getRow(v2));
        recallStat.incActualPositive();
        recallStat.incWeightActualPositive(edgeWeight);
        if (predictEdge) {
          recallStat.incTruePositive();
          recallStat.incWeightTruePositive(edgeWeight);
        }
      }
    }

    double[] colSizesNormalized = new double[matrix.getNumCols()];
    int[] colIds = new int[matrix.getNumCols()];
    double totalSize = 0;
    for (int i = 0; i < colIds.length; i++) {
      colSizesNormalized[i] = matrix.getColumn(i).size() * (matrix.getColumn(i).size() - 1);
      colIds[i] = i;
      totalSize += colSizesNormalized[i];
    }
    PredictionStat precisionStat = new PredictionStat();
    for (int i = 0; i < colIds.length; i++) {
      colSizesNormalized[i] = colSizesNormalized[i] / totalSize;
      precisionStat.add(
          clusterPrecision(
              graph, matrix, i, (int) Math.ceil(colSizesNormalized[i] * recallSamples), rng
          )
      );
    }

    PredictionStat orphansStat = new PredictionStat();
    for (int i = 0; i < verticesToSample; i++) {
      int vId = rng.nextInt(graph.getNumVertices());
      int[] row = matrix.getRow(vId);
      if (row.length > 0) {
        orphansStat.incEvalVertices();
        if (isOrphan(graph, matrix, vId, row)) {
          orphansStat.incEvalVerticesWithZeroTruePos();
        }
      }
    }

    return ImmutableMap.of(
        "precision", precisionStat,
        "recall", recallStat,
        "orphans", orphansStat
    );
  }

  /**
   * Evaluate precision of cluster i.e. what fraction of all possible edges are present
   * in the cluster. Instead of calculating by brute-force, this is calculated via sampling.
   */
  public static PredictionStat clusterPrecision(
      Graph graph,
      SparseBinaryMatrix matrix,
      int clusterId,
      int numSamples,
      RandomAdaptor rng) {
    PredictionStat ret = new PredictionStat();

    int[] colAsArray = matrix.getColumnAsArray(clusterId);
    if (colAsArray.length > 0) {
      for (int i = 0; i < numSamples;) {
        int v1 = colAsArray[rng.nextInt(colAsArray.length)];
        int v2 = colAsArray[rng.nextInt(colAsArray.length)];
        if (v1 == v2) {
          continue;
        }
        i++;
        ret.incPredictedPositive();
        float actualWeight = graph.getWeightOfEdge(v1, v2);
        if (actualWeight > 0) {
          ret.incWeightPredictedPositive(actualWeight);
          ret.incTruePositive();
          ret.incWeightTruePositive(actualWeight);
        } else {
          double add = (graph.getWeightedOutDegree(v1) + graph.getWeightedOutDegree(v2))
              / (graph.getDegree(v1) + graph.getDegree(v2));
          ret.incWeightPredictedPositive(add);
        }
      }
    }

    return ret;
  }

  /**
   * eval precision for cluster without any sampling.
   * TODO Merge with clusterPrecision above
   */
  public static void evalClusterPrecision(
      Graph graph,
      SparseBinaryMatrix matrix,
      String outputFile,
      int cpu) throws IOException {
    double[] clusterPrecision = new double[matrix.getNumCols()];
    ExecutorService executor = Executors.newFixedThreadPool(cpu);
    for (Range range : Util.chunkArray(matrix.getNumCols(), 1000)) {
      executor.submit(() -> {
        for (int k = range.start; k < range.end; k++) {
          IntSet col = matrix.getColumn(k);
          double wtTruePositive = 0;
          double wtPositive = 0;
          // SUPPRESS CHECKSTYLE NestedForDepth
          for (int v1 : col) {
            // SUPPRESS CHECKSTYLE NestedForDepth
            for (int v2 : col) {
              if (v1 < v2) {
                double edgeWt = graph.getWeightOfEdge(v1, v2);
                wtTruePositive += edgeWt;
                wtPositive += edgeWt > 0 ? edgeWt
                    : (graph.getWeightedOutDegree(v1) + graph.getWeightedOutDegree(v2))
                    / (graph.getDegree(v1) + graph.getDegree(v2));
              }
            }
          }
          clusterPrecision[k] = wtPositive > 0 ? wtTruePositive / wtPositive : 0.0;
        }
      });
    }
    executor.shutdown();
    try {
      executor.awaitTermination(60, TimeUnit.SECONDS);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
    Util.appendDoublesToFile(outputFile, clusterPrecision);
  }

  /**
   * Get precision recall numbers using vertex sampling i.e. for a subset of vertices,
   * evaluate all edges.
   */
  public static PredictionStat getPredictionStatVertexSampling(
      Graph graph,
      SparseBinaryMatrix matrix,
      ExecutorService executor,
      int[] evalVertexIds) {
    Range[] ranges = Util.chunkArray(evalVertexIds.length, 1000);
    ObjectList<Future<PredictionStat>> futures = new ObjectArrayList<>(ranges.length);
    for (Range range : ranges) {
      Future<PredictionStat> future = executor.submit(() -> {
        PredictionStat myStat = new PredictionStat();
        for (int i = range.start; i < range.end; i++) {
          int vertex1 = evalVertexIds[i];
          int truePosForVertex = 0;
          for (int vertex2 = 0; vertex2 < graph.getNumVertices(); vertex2++) {
            if (vertex1 == vertex2) {
              continue;
            }
            float edgeWeight = graph.getWeightOfEdge(vertex1, vertex2);
            boolean isEdge = edgeWeight > 0;
            boolean predictEdge =
                Util.hasCommonElement(matrix.getRow(vertex1), matrix.getRow(vertex2));
            if (isEdge) {
              myStat.incActualPositive();
              truePosForVertex++;
              myStat.incWeightActualPositive(edgeWeight);
            }
            if (predictEdge) {
              myStat.incPredictedPositive();
              // what should this be set to?
              // Let's set this to be the average weight across the two nodes incident on this edge.
              double add =
                  (graph.getWeightedOutDegree(vertex1) + graph.getWeightedOutDegree(vertex2))
                      / (graph.getDegree(vertex1) + graph.getDegree(vertex2));
              myStat.incWeightPredictedPositive(add);
            }
            if (isEdge && predictEdge) {
              myStat.incTruePositive();
              myStat.incWeightTruePositive(edgeWeight);
            }
            if (!isEdge && !predictEdge) {
              myStat.incTrueNegative();
            }
          }
          myStat.incEval(graph.getNumVertices() - 1);
          myStat.incEvalVertices();
          if (truePosForVertex == 0) {
            myStat.incEvalVerticesWithZeroTruePos();
          }
        }
        return myStat;
      });
      futures.add(future);
    }
    PredictionStat stat = new PredictionStat();
    for (Future<PredictionStat> future : futures) {
      try {
        stat.add(future.get());
      } catch (InterruptedException | ExecutionException e) {
        e.printStackTrace();
      }
    }
    return stat;
  }

  private Optional<int[]> mhStep(SparseBinaryMatrix matrix,
                                 int myId,
                                 double multiplier,
                                 double weightCoeff) {
    int[] proposal;
    if (config.proposalStrategy == ProposalStrategy.SingleOrZeroMembershipLikelihood) {
      proposal = sampleSingleNonZeroProposal(matrix, myId, weightCoeff);
    } else if (config.proposalStrategy == ProposalStrategy.MultipleMembershipLikelihood) {
      proposal =
          sampleMultipleNonZeroProposal(matrix, myId, weightCoeff, config.maxMembershipsPerVertex);
    } else if (config.proposalStrategy == ProposalStrategy.PureRandom) {
      proposal = samplePureRandomProposal(matrix.getNumCols());
    } else {
      proposal = proposalV1orV2(matrix, myId,
          config.proposalStrategy == ProposalStrategy.FractionOfNeighbors);
    }

    if (proposal.length == 0) {
      emptyProposals++;
    }
    if (!Arrays.equals(proposal, matrix.getRow(myId))) {
      // Compute likelihood
      double newLikelihood = evalLikelihood(g, matrix, myId, proposal, weightCoeff);
      double oldLikelihood = evalLikelihood(g, matrix, myId, matrix.getRow(myId), weightCoeff);

      // Accept or reject
      double logAcceptProbability = multiplier * (newLikelihood - oldLikelihood);
      if (Math.log(config.rng.nextDouble()) < logAcceptProbability) {
        return Optional.of(proposal);
      } else {
        return Optional.empty();
      }
    } else {
      sameProposals++;
      return Optional.empty();
    }
  }

  /**
   * Sum up values for each row in rowIds, from the sparse binary matrix. Keep track of
   * the non-zeros that contributed to the final result as well.
   */
  public static IdCountWtContributors[] sumRowsWithContributors(
      int[] rowIds,
      SparseBinaryMatrix matrix,
      Function<Integer, Float> rowIdToWeight) {
    ArrayList<IdCountWtContributors> ret = new ArrayList<>();

    PriorityQueue<HeapEntry> queue = new ObjectHeapPriorityQueue<>(rowIds.length);
    for (int rowId : rowIds) {
      int[] row = matrix.getRow(rowId);
      if (row.length > 0) {
        queue.enqueue(new HeapEntry(rowId, row[0], 0));
      }
    }

    double weightedSum = 0;
    int count = 0;
    IntSet contributors = new IntOpenHashSet();
    while (!queue.isEmpty()) {
      HeapEntry popped = queue.dequeue();
      weightedSum += rowIdToWeight.apply(popped.getVertexId());
      count++;
      contributors.add(popped.getVertexId());
      if (queue.isEmpty() || popped.getColId() != queue.first().getColId()) {
        ret.add(new IdCountWtContributors(popped.getColId(), count, weightedSum, contributors));
        count = 0;
        weightedSum = 0;
        contributors = new IntOpenHashSet();
      }
      int[] row = matrix.getRow(popped.getVertexId());
      if (popped.getOffset() + 1 < row.length) {
        queue.enqueue(
            new HeapEntry(
                popped.getVertexId(), row[popped.getOffset() + 1], popped.getOffset() + 1
            )
        );
      }
    }
    return Arrays.copyOf(ret.toArray(), ret.size(), IdCountWtContributors[].class);
  }

  /**
   * Sum up values for each row in rowIds, from the sparse binary matrix.
   */
  public static IdCountWtTriple[] sumRows(int[] rowIds,
                                          SparseBinaryMatrix matrix,
                                          Function<Integer, Float> rowIdToWeight) {
    ArrayList<IdCountWtTriple> ret = new ArrayList<>();

    PriorityQueue<HeapEntry> queue = new ObjectHeapPriorityQueue<>(rowIds.length);
    for (int rowId : rowIds) {
      int[] row = matrix.getRow(rowId);
      if (row.length > 0) {
        queue.enqueue(new HeapEntry(rowId, row[0], 0));
      }
    }

    double weightedSum = 0;
    int count = 0;
    while (!queue.isEmpty()) {
      HeapEntry popped = queue.dequeue();
      weightedSum += rowIdToWeight.apply(popped.getVertexId());
      count++;
      if (queue.isEmpty() || popped.getColId() != queue.first().getColId()) {
        ret.add(new IdCountWtTriple(popped.getColId(), count, weightedSum));
        count = 0;
        weightedSum = 0;
      }
      int[] row = matrix.getRow(popped.getVertexId());
      if (popped.getOffset() + 1 < row.length) {
        queue.enqueue(
            new HeapEntry(
                popped.getVertexId(), row[popped.getOffset() + 1], popped.getOffset() + 1
            )
        );
      }
    }
    return Arrays.copyOf(ret.toArray(), ret.size(), IdCountWtTriple[].class);
  }

  private int[] proposalV1orV2(
      SparseBinaryMatrix matrix,
      int vertexId,
      boolean useFractionOfNeighborsStrategy) {
    IdCountWtTriple[] factorsWithCountsWts =
        sumRows(g.getNeighbors(vertexId), matrix,
            neighborId -> g.getWeightOfEdge(vertexId, neighborId));
    int[] ret = new int[factorsWithCountsWts.length];
    int realSize = 0;
    for (IdCountWtTriple triple : factorsWithCountsWts) {
      if (triple.count <= 0) {
        throw new IllegalStateException("Triple with count zero! "
            + String.format(
            "Triple.count is %d, triple.wt is %g, triple.id is %d",
            triple.count, triple.wt, triple.id
        )
        );
      }

      double prob = 0;
      if (useFractionOfNeighborsStrategy) {
        prob = triple.wt / g.getWeightedOutDegree(vertexId);
      } else {
        double avgTruePosWeight = triple.wt / triple.count;
        double falsePosWeight =
            avgTruePosWeight * (matrix.getColumn(triple.id).size() - triple.count);
        prob = triple.wt / (g.getWeightedOutDegree(vertexId) + falsePosWeight);
      }

      if (config.rng.nextDouble() < prob) {
        ret[realSize++] = triple.id;
      }
    }

    return Arrays.copyOf(ret, realSize);
  }


  private int[] sampleMultipleNonZeroProposal(SparseBinaryMatrix matrix, int vertexId,
                                              double weightCoeff, int maxNonZero) {
    IdCountWtContributors[] sumOfRows = sumRowsWithContributors(
        g.getNeighbors(vertexId), matrix, neighborId -> g.getWeightOfEdge(vertexId, neighborId));

    double likelihoodForEmptyFactor = likelihoodInternal(g, vertexId, 0,
        0, 0, weightCoeff);

    int[] bestGumbelFactor = new int[0];
    double bestGumbelValue = likelihoodForEmptyFactor + gumbel.sample();
    double bestTruePosWt = 0;
    int bestTruePosWtId = -1;
    IntSet bestTruePostWtNeighbors = new IntOpenHashSet(0);

    for (IdCountWtContributors triple : sumOfRows) {
      if (triple.count <= 0) {
        throw new IllegalStateException("Triple with count zero! "
            + String.format(
            "Triple.count is %d, triple.wt is %g, triple.id is %d",
            triple.count, triple.wt, triple.id
        )
        );
      }

      if (config.ignoreSameProposal && matrix.getRow(vertexId).length == 1
          && matrix.getRow(vertexId)[0] == triple.id) {
        continue;
      }

      double ll = likelihoodInternal(g, vertexId, triple.count, triple.wt,
          matrix.getColumn(triple.id).size(), weightCoeff);
      if (triple.count > 0) { // only consider factors with non-zero true pos
        double gumbelValue = ll + gumbel.sample();
        if (gumbelValue > bestGumbelValue) {
          if (bestGumbelFactor.length != 1) {
            bestGumbelFactor = new int[1];
          }
          bestGumbelFactor[0] = triple.id;
          bestGumbelValue = gumbelValue;
        }
      } else {
        throw new IllegalStateException("Got triple.count zero!");
      }

      if (bestTruePosWt < triple.wt) {
        bestTruePosWt = triple.wt;
        bestTruePosWtId = triple.id;
        bestTruePostWtNeighbors = triple.contributors;
      }
    }

    if (bestTruePosWtId == -1) {
      return bestGumbelFactor;
    }

    // Keep adding columns we will consider
    int[] factor = new int[maxNonZero];
    factor[0] = bestTruePosWtId;
    IntSet factorSet = new IntOpenHashSet();
    factorSet.add(bestTruePosWtId);
    IntSet currentTruePosWtNeighbors = bestTruePostWtNeighbors;
    for (int extraIndex = 1; extraIndex < maxNonZero; extraIndex++) {
      double bestAdditionTruePosWt = 0;
      int bestAdditionTruePosWtId = -1;
      IntSet bestAdditionTruePosWtNeighbors = new IntOpenHashSet(0);

      for (IdCountWtContributors triple : sumOfRows) {
        if (!factorSet.contains(triple.id)) {
          IntSet extraContributors = new IntOpenHashSet(triple.contributors);
          extraContributors.removeAll(currentTruePosWtNeighbors);
          double wt = g.getWeightFromVertexToSet(vertexId, extraContributors);
          if (wt > bestAdditionTruePosWt) {
            bestAdditionTruePosWt = wt;
            bestAdditionTruePosWtId = triple.id;
            bestAdditionTruePosWtNeighbors = extraContributors;
          }
        }
      }

      if (bestAdditionTruePosWtId == -1) {
        break;
      }
      currentTruePosWtNeighbors.addAll(bestAdditionTruePosWtNeighbors);
      factorSet.add(bestAdditionTruePosWtId);
      factor[extraIndex] = bestAdditionTruePosWtId;
    }

    int actualMaxNonzero = factorSet.size();
    // generate all 2^actualMaxNonZero possible proposals, except for the ones we've already seen.
    for (int indicator = 0; indicator < (1 << actualMaxNonzero); indicator++) {
      if ((indicator & (indicator - 1)) == 0 || indicator == 0) {
        // let's skip powers of two, since they correspond to factors with only one factor set,
        // which we've already considered above.
        // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
        continue;
      }
      int[] newFactor = new int[maxNonZero];
      int newFactorSize = 0;
      for (int powerOf2 = 0; powerOf2 < actualMaxNonzero; powerOf2++) {
        if ((indicator & (1 << powerOf2)) > 0) {
          newFactor[newFactorSize++] = factor[powerOf2];
        }
      }

      if (config.ignoreSameProposal && Arrays.equals(matrix.getRow(vertexId), newFactor)) {
        continue;
      }

      int[] finalFactor = Arrays.copyOf(newFactor, newFactorSize);
      double ll = evalLikelihood(g, matrix, vertexId, finalFactor, weightCoeff);
      double gumbelValue = ll + gumbel.sample();
      if (gumbelValue > bestGumbelValue) {
        bestGumbelFactor = finalFactor;
        bestGumbelValue = gumbelValue;
      }
    }

    return bestGumbelFactor;
  }

  private int[] samplePureRandomProposal(int k) {
    int[] ret = new int[k];
    int realSize = 0;
    for (int i = 0; i < k; i++) {
      if (config.rng.nextBoolean()) {
        ret[realSize++] = i;
      }
    }
    return Arrays.copyOf(ret, realSize);
  }

  // http://timvieira.github.io/blog/post/2014/07/31/gumbel-max-trick/
  private int[] sampleSingleNonZeroProposal(SparseBinaryMatrix matrix, int vertexId,
                                            double weightCoeff) {
    IdCountWtTriple[] factorsWithCountsWts =
        sumRows(g.getNeighbors(vertexId), matrix,
            neighborId -> g.getWeightOfEdge(vertexId, neighborId));

    double likelihoodForEmptyFactor = likelihoodInternal(g, vertexId, 0,
        0, 0, weightCoeff);

    int bestGumbelId = -1;
    double bestGumbelValue = likelihoodForEmptyFactor + gumbel.sample();

    for (IdCountWtTriple triple : factorsWithCountsWts) {
      if (triple.count <= 0) {
        throw new IllegalStateException("Triple with count zero! "
            + String.format(
            "Triple.count is %d, triple.wt is %g, triple.id is %d",
            triple.count, triple.wt, triple.id
        )
        );
      }

      if (config.ignoreSameProposal && matrix.getRow(vertexId).length == 1
          && triple.id == matrix.getRow(vertexId)[0]) {
        // do not consider proposals which are the same as the current row for this vertex.
        continue;
      }

      double ll = likelihoodInternal(g, vertexId, triple.count, triple.wt,
          matrix.getColumn(triple.id).size(), weightCoeff);
      if (triple.count > 0) { // only consider factors with non-zero true pos
        double gumbelValue = ll + gumbel.sample();
        if (gumbelValue > bestGumbelValue) {
          bestGumbelId = triple.id;
          bestGumbelValue = gumbelValue;
        }
      } else {
        throw new IllegalStateException("Got triple.count zero!");
      }
    }

    if (bestGumbelId == -1) {
      return new int[0];
    } else {
      int[] ret = new int[1];
      ret[0] = bestGumbelId;
      return ret;
    }
  }

  private static double evalLikelihood(Graph graph, SparseBinaryMatrix matrix, int myId,
                                       int[] factor, double weightCoeff) {
    if (graph.isWeighted()) {
      return evalLikelihoodWeighted(graph, matrix, myId, factor, weightCoeff);
    } else {
      return evalLikelihoodUnweighted(graph, matrix, myId, factor, weightCoeff);
    }
  }


  private static double evalLikelihoodUnweighted(Graph graph, SparseBinaryMatrix matrix, int myId,
                                                 int[] factor, double weightCoeff) {
    int numTruePositive = 0;
    for (int neighborId : graph.getNeighbors(myId)) {
      if (Util.hasCommonElement(factor, matrix.getRow(neighborId))) {
        numTruePositive++;
      }
    }
    int numPositive = matrix.nnzInColumnSum(factor);
    int numNegative = graph.getNumVertices() - numPositive;
    int numFalseNegative = graph.getDegree(myId) - numTruePositive;
    int selfNegative = matrix.getRow(myId).length == 0 ? 1 : 0;
    int numTrueNegative = numNegative - numFalseNegative - selfNegative;
    return weightCoeff * numTruePositive + numTrueNegative;
  }

  private static double likelihoodInternal(Graph g, int vertexId, int numTruePositive,
                                           double truePosWt, int numPredictedPositive,
                                           double weightCoeff) {
    int numActualPositive = g.getDegree(vertexId);
    int numPredictedNegative = g.getNumVertices() - numPredictedPositive;
    int numFalseNegative = numActualPositive - numTruePositive;
    int numTrueNegative = numPredictedNegative - numFalseNegative;
    double trueNegativeWeight =
        numTrueNegative * g.getWeightedOutDegree(vertexId) / g.getDegree(vertexId);
    return weightCoeff * truePosWt + trueNegativeWeight;
  }

  private static double evalLikelihoodWeighted(Graph graph, SparseBinaryMatrix matrix, int vertexId,
                                               int[] factor, double weightCoeff) {
    double truePosWeight = 0;
    int numTruePositive = 0;
    if (factor.length > 0) {
      for (int neighborId : graph.getNeighbors(vertexId)) {
        if (Util.hasCommonElement(factor, matrix.getRow(neighborId))) {
          numTruePositive++;
          truePosWeight += graph.getWeightOfEdge(vertexId, neighborId);
        }
      }
    }
    return likelihoodInternal(graph, vertexId, numTruePositive, truePosWeight,
        matrix.nnzInColumnSum(factor), weightCoeff);
  }

  /**
   * This function selects a small number of neighborhoods from the input graph, based on
   * the provided cost function for a neighborhood. Note that it picks the neighborhoods
   * with the *smallest* cost function.
   */
  static IntSet[] getBottomKNeighborhoods(
      Graph g,
      int numNeighborhoods,
      BiFunction<Graph, Integer, Double> costForNeighborhoodFn,
      boolean allowOverlap,
      PrintWriter diagnosticsWriter
  ) {
    diagnosticsWriter.println("Going to get scores for all vertices");
    double tic = System.currentTimeMillis();
    int n = g.getNumVertices();
    double[] scores = new double[n];
    Integer[] vertexIds = new Integer[n];
    for (int i = 0; i < n; i++) {
      scores[i] = costForNeighborhoodFn.apply(g, i);
      vertexIds[i] = i;
      if (i > 0 && i % 1000000 == 0) {
        diagnosticsWriter.format("Done getting scores for %d vertices\r", i);
      }
    }
    // Sort scores in ascending order
    Arrays.sort(vertexIds,
        (i1, i2) -> (int) Math.signum(scores[i1] - scores[i2]));
    double toc = System.currentTimeMillis();
    diagnosticsWriter.println("Got scores for all vertices, time: " + (toc - tic) / 1000 + " secs");

    // Loop over from lowest conductance, assign each neighborhood to a community
    IntSet[] initCommunity = new IntSet[numNeighborhoods];
    IntSet seen = new IntOpenHashSet(n);
    int kSoFar = 0;
    for (int vertexId : vertexIds) {
      if (seen.contains(vertexId)) {
        continue;
      }
      IntSet c = new IntOpenHashSet(g.getDegree(vertexId));
      for (int neighbor : g.getNeighbors(vertexId)) {
        if (allowOverlap || !seen.contains(neighbor)) {
          c.add(neighbor);
        }
      }
      c.add(vertexId);
      seen.addAll(c);
      initCommunity[kSoFar++] = c;
      if (kSoFar >= numNeighborhoods) {
        break;
      }
    }
    // Assign random vertex to empty community
    while (kSoFar < numNeighborhoods) {
      IntSet randomVertexIds = new IntOpenHashSet();
      randomVertexIds.add(ThreadLocalRandom.current().nextInt(n));
      initCommunity[kSoFar++] = randomVertexIds;
    }
    return initCommunity;
  }


  /**
   * Given an list of vertex ids [0..n-1], find the index i in 0 through n - 1,
   * such that set containing vertices 0 through i has better conductance than any value
   * between 0..n-1
   * TODO move to Graph
   *
   * @param g input graph
   * @param vertexIds full set of vertexIds, ordered from most important to least important
   */
  public static int[] getBestPrefixConductanceSet(Graph g, int[] vertexIds) {
    if (vertexIds.length == 0) {
      return new int[0];
    }

    double cut = g.getWeightedOutDegree(vertexIds[0]);
    double volume = cut;
    // a set of size 1 always has worst possible conductance 1, assuming we don't count self-loops
    double bestConductanceSoFar = cut / volume;
    int argmax = 0;

    for (int i = 1; i < vertexIds.length; i++) {
      double weightedDegreeMinusCut = 0;
      for (int j = 0; j < i; j++) {
        weightedDegreeMinusCut += g.getWeightOfEdge(vertexIds[i], vertexIds[j]);
      }
      cut += g.getWeightedOutDegree(vertexIds[i]) - 2.0 * weightedDegreeMinusCut;
      volume += g.getWeightedOutDegree(vertexIds[i]);
      double newConductance =
          Math.max(cut / volume, cut / (g.getGlobalWeight() - volume));
      if (newConductance < bestConductanceSoFar) {
        bestConductanceSoFar = newConductance;
        argmax = i;
      }
    }

    return Arrays.copyOf(vertexIds, argmax + 1);
  }

  /**
   * Get best sub-neighborhoods for each vertex provided, for the given graph
   * TODO move to Graph
   */
  public static IntSet[] getBestConductanceSubNeighborhoods(Graph g, int[] vertexIds) {
    IntSet[] sets = new IntSet[vertexIds.length];
    for (int i = 0; i < vertexIds.length; i++) {
      sets[i] = new IntOpenHashSet(
          getBestPrefixConductanceSet(g, g.getNeighborsSortedByWeight(vertexIds[i]))
      );
    }

    return sets;
  }

  /**
   * Pick k random vertices and then get the best sub-neighborhoods for those vertices
   */
  public static IntSet[] getRandomBestConductanceSubNeighborhoods(Graph g, int k, Random r) {
    int[] vertexIds = new int[k];
    IntSet seen = new IntOpenHashSet();
    int n = g.getNumVertices();
    for (int i = 0; i < k; i++) {
      while (true) {
        int randomId = r.nextInt(n);
        if (!seen.contains(randomId)) {
          vertexIds[i] = randomId;
          seen.add(randomId);
          break;
        }
      }
    }
    return getBestConductanceSubNeighborhoods(g, vertexIds);
  }

  /**
   * Get k sub-neighborhoods which have good conductance and which don't overlap with each other.
   */
  public static IntSet[] getNonOverlappingBestSubNeighborhoods(Graph g, int k, Random r) {
    IntSet[] ret = new IntSet[k];
    int n = g.getNumVertices();
    IntArrayList verticesNotCovered = new IntArrayList(n);
    for (int i = 0; i < n; i++) {
      verticesNotCovered.add(i);
    }
    IntSet coveredSoFar = new IntOpenHashSet();
    for (int i = 0; i < k; i++) {
      int randomVertexId = verticesNotCovered.getInt(r.nextInt(verticesNotCovered.size()));
      int[] sortedNeighbors = g.getNeighborsSortedByWeight(randomVertexId);
      IntList unvisitedSortedNeighbors = new IntArrayList();
      for (int nId : sortedNeighbors) {
        if (!coveredSoFar.contains(nId)) {
          unvisitedSortedNeighbors.add(nId);
        }
      }
      int[] s = getBestPrefixConductanceSet(g, unvisitedSortedNeighbors.toIntArray());
      for (int j : s) {
        verticesNotCovered.rem(j);
        coveredSoFar.add(j);
      }
      ret[i] = new IntOpenHashSet(s);
    }

    return ret;
  }

  /**
   * Heuristically attach scores to the cluster assignments. The heuristic is implemented
   * in affiliationOfVertexToSet method in Graph
   * @param g graph
   * @param z sparse binary matrix with cluster assignments
   * @return Sparse real matrix with scores for each (node, cluster) pair.
   */
  public static SparseRealMatrix heuristicallyScoreClusterAssignments(Graph g,
                                                                      SparseBinaryMatrix z) {
    double[][] rowsWithScores = new double[z.getNumRows()][];
    for (int i = 0; i < z.getNumRows(); i++) {
      int[] row = z.getRow(i);
      rowsWithScores[i] = new double[row.length];
      for (int j = 0; j < row.length; j++) {
        double score = g.affiliationOfVertexToSet(i, z.getColumn(row[j]));
        rowsWithScores[i][j] = score;
      }
    }
    return new SparseRealMatrix(z, rowsWithScores);
  }
}
