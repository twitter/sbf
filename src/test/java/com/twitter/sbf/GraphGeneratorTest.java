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

package com.twitter.sbf;

import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.io.ByteStreams;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomAdaptor;
import org.junit.Test;

import com.twitter.sbf.core.AlgorithmConfig;
import com.twitter.sbf.core.MHAlgorithm;
import com.twitter.sbf.core.PredictionStat;
import com.twitter.sbf.core.ProposalStrategy;
import com.twitter.sbf.core.SparseBinaryMatrix;
import com.twitter.sbf.generator.GraphAndGroundTruth;
import com.twitter.sbf.generator.GraphGenerator;
import com.twitter.sbf.graph.Graph;

import static org.junit.Assert.assertTrue;

public class GraphGeneratorTest {

  private void checkResult(PredictionStat groundTruthResult, AlgorithmConfig cfg, Graph g,
                           SparseBinaryMatrix z, ExecutorService exec, int[] evalVertexIds) {
    checkResult(groundTruthResult, cfg, g, z, exec, evalVertexIds, g);
  }

  private void checkResult(
      PredictionStat groundTruthResult,
      AlgorithmConfig cfg,
      Graph g,
      SparseBinaryMatrix initZ,
      ExecutorService exec,
      int[] evalVertexIds,
      Graph actualGraph) {
    MHAlgorithm algo = new MHAlgorithm(cfg, g, initZ, new PrintWriter(System.err));
    SparseBinaryMatrix z = algo.optimize();
    PredictionStat algoResult =
        MHAlgorithm.getPredictionStatVertexSampling(actualGraph, z, exec, evalVertexIds);
    boolean condition = algoResult.f1() > groundTruthResult.f1() * 0.5;
    if (actualGraph.isWeighted()) {
      condition = condition && algoResult.weightedF1() > groundTruthResult.weightedF1() * 0.5;
    }
    String message = String.format("\nGround truth results:\n%s\n\n" + "Algo results:\n%s\n",
        groundTruthResult.toString(actualGraph.isWeighted()),
        algoResult.toString(actualGraph.isWeighted())
    );

    if (condition) {
      System.err.println(message);
    }
    assertTrue(message, condition);

  }

  @Test
  public void testAlgoOnWeightedGraph() {
    RandomAdaptor r = new RandomAdaptor(new JDKRandomGenerator(1));
    GraphGenerator gg = new GraphGenerator(5, 0.2,
        0.5, 0.7, 30, 60, r,
        true, 0.2f, 0.05f, 0.15f);
    GraphAndGroundTruth gAndZ = gg.generateWithDisjointClusters();
    Graph g = gAndZ.getGraph();
    SparseBinaryMatrix groundTruthZ =
        new SparseBinaryMatrix(g.getNumVertices(), gAndZ.getGroundTruth());
    System.err.println("Vertices: " + g.getNumVertices() + ", Edges: " + g.getNumEdges());
    System.err.println("Global avg. weight: " + g.getGlobalAvgWeight());

    ExecutorService exec = Executors.newFixedThreadPool(1);
    int[] evalVertexIds = g.getAllVertexIds();

    PredictionStat groundTruthResult =
        MHAlgorithm.getPredictionStatVertexSampling(g, groundTruthZ, exec, evalVertexIds);

    AlgorithmConfig cfg = new AlgorithmConfig().withK(gg.numClusters).withCpu(4).withEvalRatio(1)
        .withMaxEpoch(20).withRunAllEpochs(false).withRng(r);

    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), gg.numClusters);
    PrintWriter err = new PrintWriter(System.err);
    z.initFromBestNeighborhoods(g, (gr, i) -> cfg.rng.nextDouble(), false, err);

    checkResult(groundTruthResult, cfg, g, z, exec, evalVertexIds);
  }

  @Test
  public void testAlgoOnUnweightedGraph() {
    RandomAdaptor r = new RandomAdaptor(new JDKRandomGenerator(1));
    GraphGenerator gg = new GraphGenerator(10, 0.3,
        0.5, 0.7, 30, 80, r);
    GraphAndGroundTruth gAndZ = gg.generateWithDisjointClusters();
    Graph g = gAndZ.getGraph();
    SparseBinaryMatrix groundTruthZ =
        new SparseBinaryMatrix(g.getNumVertices(), gAndZ.getGroundTruth());
    System.err.println("Vertices: " + g.getNumVertices() + ", Edges: " + g.getNumEdges());

    ExecutorService exec = Executors.newFixedThreadPool(1);
    int[] evalVertexIds = g.getAllVertexIds();

    PredictionStat groundTruthResult =
        MHAlgorithm.getPredictionStatVertexSampling(g, groundTruthZ, exec, evalVertexIds);
    System.err.println();
    System.err.println("Ground truth result:");
    System.err.println(groundTruthResult.toString(g.isWeighted()));

    AlgorithmConfig cfg = new AlgorithmConfig().withK(gg.numClusters)
        .withEvalRatio(0.5)
        .withMaxEpoch(20)
        .withRunAllEpochs(false)
        .withRng(r);

    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), gg.numClusters);
    PrintWriter err = new PrintWriter(System.err);
    z.initFromBestNeighborhoods(
        g, (gr, i) -> cfg.rng.nextDouble(), false, err
    );
    checkResult(groundTruthResult, cfg, g, z, exec, evalVertexIds);

    exec.shutdown();
  }

  @Test
  public void testPureRandomAlgoOnGeneratedGraph() {
    RandomAdaptor r = new RandomAdaptor(new JDKRandomGenerator(1));
    GraphGenerator gg = new GraphGenerator(3, 0.3,
        0.5, 0.7, 10, 20, r);
    GraphAndGroundTruth gAndZ = gg.generateWithDisjointClusters();
    Graph g = gAndZ.getGraph();
    SparseBinaryMatrix groundTruthZ =
        new SparseBinaryMatrix(g.getNumVertices(), gAndZ.getGroundTruth());
    System.err.println("Vertices: " + g.getNumVertices() + ", Edges: " + g.getNumEdges());

    ExecutorService exec = Executors.newFixedThreadPool(1);
    int[] evalVertexIds = g.getAllVertexIds();

    PredictionStat groundTruthResult =
        MHAlgorithm.getPredictionStatVertexSampling(g, groundTruthZ, exec, evalVertexIds);
    System.err.println();
    System.err.println("Ground truth result:");
    System.err.println(groundTruthResult.toString(g.isWeighted()));

    AlgorithmConfig cfg = new AlgorithmConfig().withK(gg.numClusters)
        .withEvalRatio(0.5)
        .withMaxEpoch(20)
        .withRunAllEpochs(true)
        .withRng(r);
    cfg.proposalStrategy = ProposalStrategy.PureRandom;
    cfg.divideResultIntoConnectedComponents = false;

    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), gg.numClusters);
    z.initEmptyRowsRandomly(r);
    checkResult(groundTruthResult, cfg, g, z, exec, evalVertexIds);

    exec.shutdown();
  }

  @Test
  public void testGraphGeneration() {
    RandomAdaptor r = new RandomAdaptor(new JDKRandomGenerator(1));
    GraphGenerator gg = new GraphGenerator(10, 0.1,
        0.7, 0.9, 30, 50, r);
    System.err.format("Expected num. vertices: %s\n", gg.expectedVertices());
    System.err.format("Expected intra cluster edges: %s\n", gg.expectedIntraClusterEdges());

    GraphAndGroundTruth gAndZ = gg.generateWithDisjointClusters();
    Graph g = gAndZ.getGraph();
    System.err.println("Vertices: " + g.getNumVertices() + ", Edges: " + g.getNumEdges());

    int v = g.getNumVertices();
    long e = g.getNumEdges();

    int minVertices = gg.numClusters * gg.minClusterSize;
    int maxVertices = gg.numClusters * gg.maxClusterSize;
    assertTrue(v >= minVertices);
    assertTrue(v <= maxVertices);

    int minEdges = (int) Math.floor(
        gg.minProbInsideCluster * v * gg.minClusterSize / 2 * (1.0 + gg.fractionGlobalEdges)
    );
    int maxEdges = (int) Math.ceil(
        gg.maxProbInsideCluster * v * gg.maxClusterSize / 2 * (1.0 + gg.fractionGlobalEdges)
    );
    assertTrue(e >= minEdges);
    assertTrue(e <= maxEdges);
  }

  @Test
  public void testWeightedGraphGeneration() {
    RandomAdaptor r = new RandomAdaptor(new JDKRandomGenerator(1));
    GraphGenerator gg = new GraphGenerator(20, 0.1,
        0.7, 0.9, 30, 200, r,
        true, 0.2f, 0.05f, 0.15f);
    System.err.format("Expected num. vertices: %s\n", gg.expectedVertices());
    System.err.format("Expected intra cluster edges: %s\n", gg.expectedIntraClusterEdges());

    Graph g = gg.generateWithDisjointClusters().getGraph();
    System.err.println("Vertices: " + g.getNumVertices() + ", Edges: " + g.getNumEdges());
    System.err.println("Global avg. weight: " + g.getGlobalAvgWeight());
    System.err.println("Weighted out degree for vertex 0: " + g.getWeightedOutDegree(0)
        + ", degree: " + g.getDegree(0));

    int v = g.getNumVertices();
    long e = g.getNumEdges();

    int minVertices = gg.numClusters * gg.minClusterSize;
    int maxVertices = gg.numClusters * gg.maxClusterSize;
    assertTrue(v >= minVertices);
    assertTrue(v <= maxVertices);

    int minEdges = (int) Math.floor(
        gg.minProbInsideCluster * v * gg.minClusterSize / 2 * (1.0 + gg.fractionGlobalEdges)
    );
    int maxEdges = (int) Math.ceil(
        gg.maxProbInsideCluster * v * gg.maxClusterSize / 2 * (1.0 + gg.fractionGlobalEdges)
    );
    assertTrue(e >= minEdges);
    assertTrue(e <= maxEdges);
  }

  @Test
  public void testAlgoOnOverlappingGraph() {
    RandomAdaptor r = new RandomAdaptor(new JDKRandomGenerator(1));
    GraphGenerator gg = new GraphGenerator(20, 0.1,
        0.4, 0.7, 20, 40, r,
        true, 0.2f, 0.05f, 0.15f);

    GraphAndGroundTruth gAndZ = gg.generateWithOverlappingClusters(2,
        1.5);
    System.err.format("Expected num. vertices: %s\n", gg.expectedVertices() / 2);
    System.err.format("Expected intra cluster edges: %s\n", gg.expectedIntraClusterEdges());

    Graph g = gAndZ.getGraph();
    System.err.format("Vertices: %d, Edges: %d\n", g.getNumVertices(), g.getNumEdges());
    System.err.println("Global avg. weight: " + g.getGlobalAvgWeight());
    System.err.println("Weighted out degree for vertex 0: " + g.getWeightedOutDegree(0)
        + ", degree: " + g.getDegree(0));

    ExecutorService exec = Executors.newFixedThreadPool(1);
    int[] evalVertexIds = g.getAllVertexIds();

    SparseBinaryMatrix groundTruthZ =
        new SparseBinaryMatrix(g.getNumVertices(), gAndZ.getGroundTruth());

    PredictionStat groundTruthResult = MHAlgorithm.getPredictionStatVertexSampling(g,
        groundTruthZ, exec, evalVertexIds);

    System.err.println("Ground truth result: \n" + groundTruthResult.toString(true));
    System.err.println();

    int numClustersForAlgo = (int) Math.round(gg.numClusters * 1.5);

    AlgorithmConfig cfg = new AlgorithmConfig().withK(numClustersForAlgo).withCpu(4)
        .withEvalRatio(1)
        .withMaxEpoch(20).withRunAllEpochs(false).withRng(r);
    cfg.proposalStrategy = ProposalStrategy.MultipleMembershipLikelihood;
    cfg.maxMembershipsPerVertex = 2;

    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), numClustersForAlgo);
    z.initFromBestNeighborhoods(
        g, (gr, i) -> r.nextDouble(), false,
        new PrintWriter(new OutputStreamWriter(ByteStreams.nullOutputStream()))
    );

    System.err.println("Init from random neighborhoods: "
        + "fraction of empty rows is " + z.emptyRowProportion());
    checkResult(groundTruthResult, cfg, g, z, exec, evalVertexIds);
  }
}



