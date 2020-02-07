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

package com.twitter.sbf.simclusters_app;

import java.io.IOException;
import java.io.PrintWriter;


import com.twitter.sbf.core.AlgorithmConfig;
import com.twitter.sbf.core.BipartiteSparseRepresentations;
import com.twitter.sbf.core.SimClustersInMemory;
import com.twitter.sbf.core.SparseRealMatrix;
import com.twitter.sbf.graph.BipartiteGraph;

public final class Main {
  private Main() {
  }

  private static ConfigSimClusters config;

  /**
   * Run as stand-alone
   */
  public static void main(String[] args) throws IOException {
    //If config file doesn't exist, bust
    if (args.length != 1) {
      throw new IllegalStateException("Usage: command config_filename");
    }
    //Load config data into object
    config = new ConfigSimClusters(args[0]);
    //Run SimClusters
    run();
  }

  private static void run() throws IOException {

    // Load graph from file
    System.out.println("Loading bipartite graph from METIS data");
    BipartiteGraph graph = BipartiteGraph.fromFile(config.getMetisFile());
    System.out.println("Load graph: done!");

    //Variable for specifying the type of initialization for the algorithm
    final String evalType;
    if (config.getInitFromRandomNeighborhoods()) {
      evalType = "randomNeighborhoods";
    } else if (config.getInitFromBestNeighborhood()) {
      evalType = "bestNeighborhoods";
    } else if (config.getInitFromNonoverlappingNeighborhood()) {
      evalType = "nonoverlappingNeighborhoods";
    } else {
      evalType = "random";
    }

    //Arguments for SimClusters object
    int k = config.getK();
    AlgorithmConfig paramsForStepTwo = config.getAlgoConfig();
    int thresholdForStepThree = config.getThresholdForStepThree();
    double thresholdForStepOne = config.getThresholdForStepOne();
    boolean retainTopK = config.getRetainTopK();
    int topKParameter = config.getTopKParameter();
    boolean normalizeColumnsRight = config.getNormalizeColumnsRight();
    boolean normalizeRowsLeft = config.getNormalizeRowsLeft();
    boolean applyLogTransform = config.getApplyLogTransform();
    boolean squareWeights = config.getSquareWeights();
    boolean applyStepFour = config.getApplyStepFour();
    double thresholdForStepFour = config.getThresholdForStepFour();
    PrintWriter log = new PrintWriter(System.out);
    //Create SimClusters algorithm object
    SimClustersInMemory algo
        = new SimClustersInMemory(
            graph, k, paramsForStepTwo, thresholdForStepThree, thresholdForStepOne, retainTopK,
        topKParameter, log, normalizeColumnsRight, normalizeRowsLeft, squareWeights, evalType,
        applyStepFour, thresholdForStepFour, applyLogTransform);
    //Run algorithm
    long tic = System.currentTimeMillis();
    BipartiteSparseRepresentations representations =
        algo.runSimClusters();
    long toc = System.currentTimeMillis();

    System.out.println(String.format(
        "Time to learn bipartite sparse representations: %.2f seconds", (toc - tic) / 1000.0
        ));

    // Write output
    System.out.println(
        "Printing left representations to file " + config.getOutputLeftRepFile()
        + " and right representations to file " + config.getOutputRightRepFile()
    );
    PrintWriter leftOutput = new PrintWriter(config.getOutputLeftRepFile());
    PrintWriter rightOutput = new PrintWriter(config.getOutputRightRepFile());
    writeRowsWithScores(representations.leftRepresentations, leftOutput);
    writeRowsWithScores(representations.rightRepresentations, rightOutput);
    System.out.println("Done writing output!");
  }

  /**
   * Function for writing a sparse real matrix to a file. The format of the output is as follows:
   * - The first line contains "m n" (the dimensions of the matrix)
   * - The next m lines describe each of the m rows specifying the column ids for the non-zero
   *   columns along with the double value associated with the entry
   *
   * @param z - SparseRealMatrix to be printed to file
   * @param writer - A PrintWriter that writes to the file
   * @throws IOException
   */
  private static void writeRowsWithScores(
      SparseRealMatrix z,
      PrintWriter writer
  ) throws IOException {
    //Print first line
    writer.println(String.format("%d %d", z.getNumRows(), z.getNumCols()));
    //Print row-wise data
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



