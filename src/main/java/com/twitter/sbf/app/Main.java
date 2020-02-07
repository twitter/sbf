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
 
package com.twitter.sbf.app;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import com.twitter.sbf.core.MHAlgorithm;
import com.twitter.sbf.core.PredictionStat;
import com.twitter.sbf.core.SparseBinaryMatrix;
import com.twitter.sbf.core.SparseRealMatrix;
import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.util.Util;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

public final class Main {
  private Main() {
  }

  private static Config config;

  /**
   * Run as stand-alone
   */
  public static void main(String[] args) throws IOException {
    if (args.length != 1) {
      System.err.println("Usage: command config_filename");
    } else {
      config = new Config(args[0]);
      System.err.println(config.toString());

      run();
    }
  }

  private static void run() throws IOException {
    // Sanity check
    if (!config.clusterPrecisionFile.isEmpty() && Util.fileExists(config.clusterPrecisionFile)) {
      throw new IllegalStateException("clusterPrecisionFile already exists!");
    }

    // Load graph from file
    System.out.println("Loading graph");
    Graph graph = Graph.fromFile(config.metisFile);
    System.out.println("Load graph: done");
    graph.print();

    // Allocate empty Z
    int n = graph.getNumVertices();
    SparseBinaryMatrix z = new SparseBinaryMatrix(n, config.getAlgoConfig().k);
    PrintWriter err = new PrintWriter(System.err);

    long tic = System.currentTimeMillis();
    // Initialize Z
    if (!config.initFromRowsFile.isEmpty()) {
      System.out.println("Initializing factors from rows");
      z.initFromRows(config.initFromRowsFile);
      System.out.println("Initialization from rows: done");
    } else if (!config.initFromColsFile.isEmpty()) {
      System.out.println("Initializing factors from columns");
      IntSet[] initCols = readColumnsFromFile(config.initFromColsFile);
      if (initCols.length != config.getAlgoConfig().k) {
        System.out.format(
            "Number of columns in %d different from K specified in config %d\n",
            initCols.length, config.getAlgoConfig().k
        );
        System.out.println("Will use the number of columns in the file as the new K");
        z = new SparseBinaryMatrix(n, initCols.length);
      }
      z.initFromColSets(initCols);
      System.out.println("Initialization from columns: done");
    } else if (config.initFromRandomNeighborhoods) {
      System.out.println("Initializing from random neighborhoods");
      // set allowOverlap = false
      z.initFromBestNeighborhoods(
          graph, (g, i) -> config.getAlgoConfig().rng.nextDouble(), false, err
      );
      PredictionStat prec0 =
          MHAlgorithm.clusterPrecision(
              graph, z, 0, 1000, config.getAlgoConfig().rng
          );
      System.out.println("Precision of cluster 0:" + prec0.precision());
      PredictionStat prec1 =
          MHAlgorithm.clusterPrecision(
              graph, z, 1, 1000, config.getAlgoConfig().rng
          );
      System.out.println("Precision of cluster 1:" + prec1.precision());
      System.out.println(
          "Fraction of empty rows after initializing from random neighborhoods: "
              + z.emptyRowProportion()
      );
    } else if (config.initFromBestNeighborhood) {
      System.out.println("Initializing from best sub-neighborhoods in terms of conductance");
      z.initFromColSets(
          MHAlgorithm.getRandomBestConductanceSubNeighborhoods(
              graph, z.getNumCols(), config.getAlgoConfig().rng
          )
      );
      System.out.println(
          "Fraction of empty rows after initializing using best sub-neighborhoods: "
              + z.emptyRowProportion()
      );
      //Z.initEmptyRowsRandomly(config.getAlgoConfig().rng);
      //System.out.println("Initializing from best neighborhoods in terms of conductance: done");
    } else if (config.initFromNonoverlappingNeighborhood) {
      System.out.println(
          "Initializing from best non-overlapping sub-neighborhoods in terms of conductance");
      z.initFromColSets(
          MHAlgorithm.getNonOverlappingBestSubNeighborhoods(
              graph, z.getNumCols(), config.getAlgoConfig().rng
          )
      );
      System.out.println(
          "Fraction of empty rows after initializing using best non-overlapping sub-neighborhoods: "
              + z.emptyRowProportion()
      );
    } else {
      System.out.println("Initializing factors randomly with 1 nnz/vertex");
      z.initEmptyRowsRandomly(config.getAlgoConfig().rng);
      System.out.println("Random initialization: done");
    }

    long toc = System.currentTimeMillis();
    System.out.println(String.format("Time to initialize: %.2f seconds\n", (toc - tic) / 1000.0));
    MHAlgorithm algo = new MHAlgorithm(config.getAlgoConfig(), graph, z, err);
    SparseBinaryMatrix optimizedZ = algo.optimize();
    long toc2 = System.currentTimeMillis();
    System.out.println(String.format("Time to optimize: %.2f seconds\n", (toc2 - toc) / 1000.0));
    System.out.println(String.format("Time to initialize & optimize: %.2f seconds\n",
      (toc2 - tic) / 1000.0));

    // Write output
    if (!config.outputByRowsFile.isEmpty()) {
      System.out.println("Writing rows of Z to " + config.outputByRowsFile);
      optimizedZ.outputByRows(config.outputByRowsFile);
      System.out.println("Output by rows: done");
    }
    if (!config.outputByColsFile.isEmpty()) {
      System.out.println("Writing columns of Z to " + config.outputByColsFile);
      optimizedZ.outputByCols(config.outputByColsFile);
      System.out.println("Output by columns: done");
    }
    if (!config.clusterPrecisionFile.isEmpty()) {
      MHAlgorithm.evalClusterPrecision(graph, optimizedZ,
          config.clusterPrecisionFile, config.getAlgoConfig().cpu);
    }
    if (!config.outputRowsWithScoresFile.isEmpty()) {
      System.out.println("Writing rows of Z with scores to " + config.outputRowsWithScoresFile);
      SparseRealMatrix srm =
          MHAlgorithm.heuristicallyScoreClusterAssignments(graph, optimizedZ);
      srm.normalizeToUnitColumn();
      PrintWriter w = new PrintWriter(config.outputRowsWithScoresFile);
      writeRowsWithScores(srm, w);
      System.out.println("Output rows with scores: done");
    }
  }

  /**
   * Read from file that contains column view of the matrix.
   * Each line contains nonzero row id's of a column, separated by whitespace.
   * NOTE: assume row id's are 1-based, to accommodate output from SNAP.
   *
   * @param filename (required) file where each line is a column
   */
  private static IntSet[] readColumnsFromFile(String filename) throws IOException {
    int numLines = (int) Util.getNumLinesOfFile(filename);
    IntSet[] columns = new IntSet[numLines];
    // Read from file
    BufferedReader br = new BufferedReader(new FileReader(filename));
    int colId = 0;
    String line;
    while ((line = br.readLine()) != null) {
      line = line.trim();
      if (!line.isEmpty()) {
        String[] tokens = line.split("\\s+");
        columns[colId] = new IntOpenHashSet(tokens.length);
        for (String token : tokens) {
          int rowId = Integer.parseInt(token) - 1; // row id's are 1-based, so make offset
          columns[colId].add(rowId);
        }
      } else {
        columns[colId] = new IntOpenHashSet(0);
      }
      colId++;
    }
    br.close();
    return columns;
  }

  private static void writeRowsWithScores(
      SparseRealMatrix z,
      PrintWriter writer
  ) throws IOException {
    for (int i = 0; i < z.getNumRows(); i++) {
      int[] rowWithIndices = z.getColIdsForRow(i);
      double[] rowWithScores = z.getValuesForRow(i);
      for (int j = 0; j < rowWithIndices.length; j++) {
        if (j > 0) {
          writer.print(" ");
        }
        // Add 1 to get back 1-indexing
        writer.print(String.format("%d:%.2g", rowWithIndices[j] + 1, rowWithScores[j]));
      }
      writer.println();
    }
    writer.close();
  }

}
