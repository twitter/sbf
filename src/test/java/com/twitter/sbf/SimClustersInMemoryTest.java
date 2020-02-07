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

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomAdaptor;
import org.junit.Test;

import com.twitter.sbf.core.AlgorithmConfig;
import com.twitter.sbf.core.BipartiteSparseRepresentations;
import com.twitter.sbf.core.MHAlgorithm;
import com.twitter.sbf.core.SimClustersInMemory;
import com.twitter.sbf.core.SparseBinaryMatrix;
import com.twitter.sbf.core.SparseRealMatrix;
import com.twitter.sbf.graph.BipartiteGraph;
import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.util.SimpleIteratorFromIterator;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

public class SimClustersInMemoryTest {
  //Returns bipartite graph with 15 vertices on left, 9 vertices on the right, and 42 edges in total
  //in the METIS format.
  //Graph has 3 bipartite communities
  private static BipartiteGraph loadEasyBipartiteGraph() {
    List<String> lines = new ArrayList<>();
    lines.add("15 9 42");
    //Adjacency lists for left side vertices
    lines.add("1 2 3"); // 1
    lines.add("1 2 3"); // 2
    lines.add("1 2 3"); // 3
    lines.add("1 2 3"); // 4
    lines.add("1 4"); // 5
    lines.add("4 5 6"); // 6
    lines.add("4 5 6"); // 7
    lines.add("4 5 6"); // 8
    lines.add("4 5 6"); // 9
    lines.add("4 7"); // 10
    lines.add("7 8 9"); // 11
    lines.add("7 8 9"); // 12
    lines.add("7 8 9"); // 13
    lines.add("7 8 9"); // 14
    lines.add("1 7"); // 15
    //Adjacency lists for right side vertices
    lines.add("1 2 3 4 5 15"); // 1
    lines.add("1 2 3 4"); // 2
    lines.add("1 2 3 4"); // 3
    lines.add("5 6 7 8 9 10"); // 4
    lines.add("6 7 8 9"); // 5
    lines.add("6 7 8 9"); //6
    lines.add("10 11 12 13 14 15"); // 7
    lines.add("11 12 13 14"); // 8
    lines.add("11 12 13 14"); // 9
    return new BipartiteGraph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  //Returns an AlgorithmConfig object pre-configured for unit tests
  private static AlgorithmConfig loadConfigForMH() {
    AlgorithmConfig cfg =
        new AlgorithmConfig().withK(3).withScaleCoeff(1.0).withWtCoeff(5.0)
            .withRunAllEpochs(true).withMaxEpoch(10).withCpu(1).withEvalRatio(1.0)
            .withRng(new RandomAdaptor(new JDKRandomGenerator(1)));
    cfg.updateImmediately = true;
    return cfg;
  }

  /**
   * This method returns a SparseBinaryMatrix objects that returns an initialized matrix to be used
   * for the Metropolis-Hastings step of SimClusters. Uses neighbourhood conductance for
   * initialization.
   * Returns a initialized SparseBinaryMatrix object of dimensions numVertices in g X k
   *
   * @param g The undirected graph which is to be used for initialization
   * @param k The number of communities/latent factors
   * @return SparseBinaryMatrix object
   */
  private static SparseBinaryMatrix loadInitializedMatrix(Graph g, int k) {
    SparseBinaryMatrix sIn = new SparseBinaryMatrix(g.getNumVertices(), k);
    PrintWriter err = new PrintWriter(System.err);
    sIn.initFromBestNeighborhoods(g, Graph::getNeighborhoodConductance, false, err);
    return sIn;
  }

  /**
   * Method that returns a randomly initialized SparseBinaryMatrix for Metropolis-Hastings step
   * of SimClusters
   *
   * @param g A Graph object
   * @param cfg An AlgorithmConfig object
   * @return SparseBinaryMatrix object
   */
  private static SparseBinaryMatrix loadRandomlyInitializedMatrix(Graph g, AlgorithmConfig cfg) {
    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), cfg.k);
    z.initEmptyRowsRandomly(cfg.rng);
    return z;
  }

  /**
   * Method to print a Graph object (weighted)
   */
  private static void printGraph(Graph g) {
    int n = g.getNumVertices();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        System.out.print(String.format("%f ", g.getWeightOfEdge(i, j)));
      }
      System.out.println(" ");
    }
  }

  /**
   * Method to print a SparseRealMatrix object
   *
   * @param matrix - SparseRealMatrix object
   */
  private static void printSparseRealMatrix(SparseRealMatrix matrix) {
    for (int i = 0; i < matrix.getNumRows(); i++) {
      int[] nzInCurrentRow = matrix.getColIdsForRow(i);
      double[] valuesInCurrentRow = matrix.getValuesForRow(i);
      int head = 0;
      for (int curCol = 0; curCol < matrix.getNumCols(); curCol++) {
        if (nzInCurrentRow.length > 0 && curCol == nzInCurrentRow[head]) {
          System.out.print(String.format("%f ", valuesInCurrentRow[head]));
          head = (head < nzInCurrentRow.length - 1) ? head + 1 : nzInCurrentRow.length - 1;
        } else {
          System.out.print("0 ");
        }
      }
      System.out.print("\n");
    }
  }

  /**
   * Method to check if two given SparseRealMatrix objects are equal
   */
  private void testMatrixEquality(SparseRealMatrix matrix1, SparseRealMatrix matrix2) {
    //Check if dimensions match
    assertTrue(matrix1.getNumCols() == matrix2.getNumCols());
    assertTrue(matrix1.getNumRows() == matrix2.getNumRows());

    //Check if the support of each row is the same in both matrices
    for (int rowId = 0; rowId < matrix1.getNumRows(); rowId++) {
      assertArrayEquals(matrix1.getColIdsForRow(rowId), matrix2.getColIdsForRow(rowId));
    }

    //Finally check if the non-zero values are same in each row
    for (int rowId = 0; rowId < matrix1.getNumRows(); rowId++) {
      double[] nzInCurRow1 = matrix1.getValuesForRow(rowId);
      double[] nzInCurRow2 = matrix2.getValuesForRow(rowId);
      for (int colId = 0; colId < nzInCurRow1.length; colId++) {
        assertTrue(Math.abs(nzInCurRow1[colId] - nzInCurRow2[colId]) < 1e-5);
      }
    }
  }

  //Method to load the expected representations for the right vertices
  //We use threshold = 0.0 for the project step. See BipartiteGraph for more details.
  private static SparseRealMatrix getExpectedRightRepresentations() {
    //Load example graph
    BipartiteGraph g = loadEasyBipartiteGraph();
    //project graph on the right side
    Graph projectedGraph = g.projectRight(0.0);
    //load configuration
    AlgorithmConfig cfg = loadConfigForMH();
    //load initialized matrix
    SparseBinaryMatrix z = loadInitializedMatrix(projectedGraph, 3);
    MHAlgorithm algo = new MHAlgorithm(cfg, projectedGraph, z, new PrintWriter(System.err));
    //Run MHAlgorithm
    SparseBinaryMatrix s = algo.optimize();
    SparseRealMatrix expectedOutput =
        MHAlgorithm.heuristicallyScoreClusterAssignments(projectedGraph, s);
    //expectedOutput is the output we expect when we run getRepresentationsForRight().
    //The step below has been commented out. If you need every column on the representationsForRight
    //matrix to be normalized to unit then please uncomment it.
    //expectedOutput.normalizeToUnitColumn();
    //return expected matrix
    return expectedOutput;
  }

  //Method to load the expected representations for the left vertices
  private static SparseRealMatrix getExpectedLeftRepresentations() {
    //Specify the 3 columns as IntSet[]
    int[] col1 = {9, 10, 11, 12, 13, 14};
    int[] col2 = {0, 1, 2, 3, 4, 14};
    int[] col3 = {4, 5, 6, 7, 8, 9};
    int[][] colArrays = new int[3][];
    colArrays[0] = col1;
    colArrays[1] = col2;
    colArrays[2] = col3;
    //Convert above columns to IntSets
    IntSet[] cols = new IntSet[3];
    for (int i = 0; i < 3; i++) {
      cols[i] = new IntOpenHashSet();
      for (int index : colArrays[i]) {
        cols[i].add(index);
      }
    }
    SparseBinaryMatrix expectedSupport = new SparseBinaryMatrix(15, 3);
    expectedSupport.initFromColSets(cols);
    //Specify the non-zero values in the rows
    //There are only two distinct non-zero values that occur in the matrix
    double valBig = 2.830479;
    double valSmall = 0.830479;
    double[][] rowVals = new double[15][];
    //rows that only have one non-zero entry equal to valBig
    int[] singleSupportRowIndices = {0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13};
    //rows that have two non-zero entries both equal to valSmall
    int[] doubleSupportRowIndices = {4, 9, 14};
    //Initialize the single support rows first
    for (int index : singleSupportRowIndices) {
      rowVals[index] = new double[1];
      rowVals[index][0] = valBig;
    }
    //Initialize the double support rows next
    for (int index : doubleSupportRowIndices) {
      rowVals[index] = new double[2];
      rowVals[index][0] = valSmall;
      rowVals[index][1] = valSmall;
    }
    //define and return SparseRealMatrix
    return new SparseRealMatrix(expectedSupport, rowVals);
  }

  //Method to load the expected representations for the left vertices when thresholdForStepThree = 2
  private static SparseRealMatrix getExpectedLeftRepresentationsWithThreshold() {
    //Specify the 3 columns as IntSet[]
    int[] col1 = {10, 11, 12, 13};
    int[] col2 = {0, 1, 2, 3};
    int[] col3 = {5, 6, 7, 8};
    int[][] colArrays = new int[3][];
    colArrays[0] = col1;
    colArrays[1] = col2;
    colArrays[2] = col3;
    //Convert above columns to IntSets
    IntSet[] cols = new IntSet[3];
    for (int i = 0; i < 3; i++) {
      cols[i] = new IntOpenHashSet();
      for (int index : colArrays[i]) {
        cols[i].add(index);
      }
    }
    SparseBinaryMatrix expectedSupport = new SparseBinaryMatrix(15, 3);
    expectedSupport.initFromColSets(cols);
    //Specify the non-zero values in the rows
    //There are only two distinct non-zero values that occur in the matrix
    double valBig = 2.830479;
    double[][] rowVals = new double[15][];
    //rows that only have one non-zero entry equal to valBig
    int[] singleSupportRowIndices = {0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13};
    //Initialize the single support rows first
    for (int index : singleSupportRowIndices) {
      rowVals[index] = new double[1];
      rowVals[index][0] = valBig;
    }
    //define and return SparseRealMatrix
    return new SparseRealMatrix(expectedSupport, rowVals);
  }

  //Method to load expected right representations matrix after applying step 4 without thresholding
  private static SparseRealMatrix getExpectedUpdatedRightRepresentations() {
    //Specify support of 3 columns as IntSet[]
    int[] col1 = {0, 3, 6, 7, 8};
    int[] col2 = {0, 1, 2, 3, 6};
    int[] col3 = {0, 3, 4, 5, 6};
    int[][] colArrays = new int[3][];
    colArrays[0] = col1;
    colArrays[1] = col2;
    colArrays[2] = col3;
    //Convert above columns to IntSets
    IntSet[] cols = new IntSet[3];
    for (int i = 0; i < 3; i++) {
      cols[i] = new IntOpenHashSet();
      for (int index : colArrays[i]) {
        cols[i].add(index);
      }
    }
    //Create non-zero support structure in the form of SparseBinaryMatrix
    SparseBinaryMatrix expectedSupport = new SparseBinaryMatrix(9, 3);
    expectedSupport.initFromColSets(cols);
    //Specify values in rows as double[][]
    double[][] rowVals = new double[9][];
    double[] row1 = {0.058642420729399944012, 0.91675666621377815982, 0.058642420729399944012};
    double[] row2 = {0.9791490172};
    double[] row3 = {0.9791490172};
    double[] row4 = {0.058642420729399944012, 0.058642420729399944012, 0.91675666621377815982};
    double[] row5 = {0.9791490172};
    double[] row6 = {0.9791490172};
    double[] row7 = {0.91675666621377815982, 0.058642420729399944012, 0.058642420729399944012};
    double[] row8 = {0.9791490172};
    double[] row9 = {0.9791490172};
    rowVals[0] = row1;
    rowVals[1] = row2;
    rowVals[2] = row3;
    rowVals[3] = row4;
    rowVals[4] = row5;
    rowVals[5] = row6;
    rowVals[6] = row7;
    rowVals[7] = row8;
    rowVals[8] = row9;
    //create and return new SparseRealMatrix based on above
    return new SparseRealMatrix(expectedSupport, rowVals);
  }

  @Test
  //Method to test the runSimClusters method of the class
  public void testRunSimClusters() {
    //Load example bipartite graph
    BipartiteGraph g = loadEasyBipartiteGraph();
    //Get projection for initialization
    Graph projectedGraph = g.projectRight(0.0);
    //Get the default configuration and initialization matrix
    AlgorithmConfig testConfig = loadConfigForMH();
    //Define new SimClustersInMemory object
    SimClustersInMemory testAlg =
        new SimClustersInMemory(
            g, 3, testConfig, 1, 0.0, false,
            -1, new PrintWriter(System.err), false,
            false, false, "test",
            false, 0.0, false
        );
    //Run SimClusters to obtain representations for left and right vertices.
    BipartiteSparseRepresentations outputOfSimClusters = testAlg.runSimClusters();
    //Test if left and right representations are as expected
    testMatrixEquality(outputOfSimClusters.leftRepresentations, getExpectedLeftRepresentations());
    testMatrixEquality(outputOfSimClusters.rightRepresentations, getExpectedRightRepresentations());

  }

  @Test
  //Method to test the runSimClusters method with thresholdForStepThree = 2
  public void testRunSimClustersWithThresholding() {
    //Load example bipartite graph
    BipartiteGraph g = loadEasyBipartiteGraph();
    //Get projection for initialization
    Graph projectedGraph = g.projectRight(0.0);
    //Get the default configuration and initialization matrix
    AlgorithmConfig testConfig = loadConfigForMH();
    //Repeat tests with thresholdForStepThree = 2 keeping everything else the same
    SimClustersInMemory testAlgWithThreshold2 =
        new SimClustersInMemory(
            g, 3, testConfig, 2, 0.0, false,
            -1, new PrintWriter(System.err), false,
            false, false, "test",
            false, 0.0, false
        );
    BipartiteSparseRepresentations outputOfSimClustersWithThreshold =
        testAlgWithThreshold2.runSimClusters();
    //There should be no change in the right representations
    testMatrixEquality(
        outputOfSimClustersWithThreshold.rightRepresentations, getExpectedRightRepresentations()
    );
    //As for left representations, vertex 5, 10, and 15 (in 1-based indexing) on the left do not
    //any non-zero value in their latent representation.
    testMatrixEquality(
        outputOfSimClustersWithThreshold.leftRepresentations,
        getExpectedLeftRepresentationsWithThreshold()
    );
  }

  @Test
  //Method to test SimClusters with applyStepFour = true and thresholdForStepFour =
  public void testRunSimClustersWithStepFour() {
    //Load example bipartite graph
    BipartiteGraph g = loadEasyBipartiteGraph();
    //The expected right representation matrix after applying step four
    SparseRealMatrix expectedUpdatedRightRepresentations = getExpectedUpdatedRightRepresentations();
    //Get default configuration for step two
    AlgorithmConfig testConfig = loadConfigForMH();
    //Create new SimClustersInMemory object with the all the right parameters
    SimClustersInMemory testAlgWithStepFour =
        new SimClustersInMemory(
            g, 3, testConfig, 0, 0.0, false,
            -1, new PrintWriter(System.err), false,
            false, false, "test",
            true, 0.0, false
        );
    //Run algorithm and get representations for left and right
    BipartiteSparseRepresentations outputOfSimClustersWithStepFour =
        testAlgWithStepFour.runSimClusters();
    printSparseRealMatrix(outputOfSimClustersWithStepFour.rightRepresentations);
    //Left representations have already been compared and tested in other tests; check right
    testMatrixEquality(
        outputOfSimClustersWithStepFour.rightRepresentations,
        expectedUpdatedRightRepresentations
    );
  }
}
