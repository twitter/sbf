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

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Set;

import com.twitter.sbf.graph.BipartiteGraph;
import com.twitter.sbf.graph.Graph;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;


/**
 * Class that performs implements the SimClusters in-memory algorithm. Given an input
 * bipartite graph g with m left vertices and n right vertices, and an input parameter k equal
 * to a guess for the number of communities (number of latent factors), the algorithm returns m+n
 * non-negative vectors of dimension k' (Note k' could be different from k!)
 * A representation vector of dimension k' of a vertex can be thought of as the strengths of
 * affiliation of the vertex to the k' communities. The algorithm proceeds in 3 stages and assumes
 * that the caller provides the class with an initialization matrix. The algorithm has an optional
 * 4th step that can be performed if required and specified by the caller.
 * <p>
 * 1.) Without loss of generality assume that m > n. The algorithm performs a unipartite weighted
 * projection on the input bipartite as in the method BipartiteGraph.projectRight(). This
 * results in a weighted undirected graph on n vertices where the edge weights capture
 * similarity between pairs of right side vertices (similarity as far as their neighborhoods in
 * the bipartite graphs are concerned). We also perform thresholding during this step to prune some
 * edges. See BipartiteGraph.projectRight() for more details.
 * There is also an additional option of retaining only the top-k edges for every vertex after
 * projection and of squaring the weights. (See the class Graph for more details)
 * (This is a different k from the dimension)
 * <p>
 * 2.) The algorithm then performs community detection using Metropolis-Hastings sampling on the
 * weighted undirected graph obtained in step 1. This is implemented in the class MHAlgorithm.
 * The result is vector of dimension k' (k' would be different from k) for each of the n
 * vertices representing the strengths of affiliation of a vertex to the k' communities. The
 * design of the Metropolis-Hastings based algorithm is such that the vector for every vertex
 * has at most one non-negative entry, i.e. every vertex on the right is affiliated with at
 * most one community. Note again that k' is not necessarily equal to k.
 * <p>
 * 3.) Let us denote the m X n adjacency matrix of the input bipartite graph by A and that of the
 * weighted undirected graph on n vertices obtain at the end of step 1 by S. Step 2 can
 * essentially be thought of as approx. finding a non-negative matrix of dimension n X k', say V,
 * such that S = VV^T, where V^T is the transpose of V.
 * In Step 3, the algorithm basically computes AV to obtain a non-negative m X k' matrix that
 * gives us the representation vectors for the m vertices on the left.
 * The algorithm returns two SparseRealMatrices bundled together as a BipartiteSparseRepresentations
 * object.
 * Additionally, we also implement the following thresholding operation at this step:
 * - Let 1 <= i <= m be the index of a vertex on the left hand side and let
 * 1 <= j < = k' be the index of a latent feature/community. Then we set (i,j) entry of the product
 * matrix AV to be zero if vertex i is connected to at most one right vertex that belongs to
 * community j (i.e., has a non-zero value in index j in its latent representation).
 * There is also the option of applying a coordinate-wise log10(1+x) tranform at the end of Step 3.
 * 4.) (Optional step***) There is a variable "applyStepFour", which if set to True, specifies that
 * the fourth step of the algorithm also needs to be performed. In this step, we update the
 * representation vectors for the left matrices as follows:
 * If U is the m X k' left representations matrix, and A is the m X n bi-adjacency matrix, then
 * in the updated right representations matrix V_updated, row i and column j contains the value:
 * cosine(column i of A, column j of U).
 * Thus, V_updated is an n X k' matrix.
 * <p>
 * NOTE: We expect all the graphs and matrices involved in the algorithm to be sparse (especially
 * for large m and n) and so we use sparse data structures and algorithms for sparse matrices
 * and data structures to implement the algorithm.
 * <p>
 * NOTE: This implementation only works for unweighted bipartite graphs. If in the future the
 * BipartiteGraph class is extended to weighted bipartite graphs, then the code below would need to
 * be modified.
 */
public class SimClustersInMemory {
  //The input bipartite graph of dimensions m X n, m left and n right vertices.
  private BipartiteGraph g;
  //Guess for number of communities or latent factors per vertex that the caller expects there to be
  //in the graph.
  private int k;
  //An AlgorithmConfig that stores the parameters for Step 2 of the algorithm, i.e.
  //the Metropolis-Hastings algorithm. See MHAlgorithm.java and AlgorithmConfig.java for details.
  private AlgorithmConfig mhParameters;
  //The thresholding parameter for step three of the algorithm (See class description given above)
  private int thresholdForStepThree;
  //The thresholding parameter for step one (BigraphGraph.projectRight(). See function for details)
  private double thresholdForStepOne;
  //Boolean parameter to specify whether to only keep top-k edges for every vertex after projection.
  //This is a different k from the dimension k above
  private boolean retainTopK;
  //If retainTopK is true, then what is the value of k
  private int topKParameter;
  //A PrintWriter object for diagnostic messages.
  private PrintWriter diagnosticsWriter;
  //Should the columns of right representation matrix be normalized
  private boolean normalizeColumnsRight;
  //Should the rows of the left representations matrix be normalized
  private boolean normalizeRowsLeft;
  //Variable for specifying if weights should be squared in first step.
  private boolean squareWeights;
  //Variable for specifying type of init. for second step.,
  private String initType;
  //Boolean for specifying if step four needs to be performed
  private boolean  applyStepFour;
  //The thresholding parameter for step four of the algorithm (See class description) if it is
  //to be performed
  private double thresholdForStepFour;
  //Boolean for specifying if log10(1+x) coordinate-wise tranform should be applied after Step 3
  private boolean applyLogTransform;

  /**
   * Constructor method
   * This method sets defaults value for the AlgorithmConfig object except for the parameter k
   *
   * @param g - a BipartiteGraph object representing the graph
   * @param k - the number of communities that the caller expects there to be
   * @param mhParameters - parameters for the step 2 of the algorithm, i.e. Metropolis-Hastings
   * @param thresholdForStepThree - thresholding parameter for step 3 (see above)
   * @param thresholdForStepOne - thresholding parameter for step 1 (see above)
   * @param retainTopK - whether to perform top-k operation after projection
   * @param topKParameter - if top-k retention is enabled then what is k
   * @param diagnosticsWriter - a PrintWriter object for writing diagnostics messages
   * @param normalizeColumnsRight - whether columns should be normalized in right matrix
   * @param normalizeRowsLeft - whether rows in left matrix should be normalized
   * @param squareWeights - whether weights should be squared in projected graph.
   * @param initType - what type of init. to use for step 2
   * @param applyStepFour - whether or not to perform step four of the algorithm.
   * @param thresholdForStepFour - the threshold value for step 4 (see above).
   * @param applyLogTransform - whether to apply log10(1+x) transform after Step 3
   * NOTE: The k in the AlgorithmConfig object must match the k passed as argument.
   */
  public SimClustersInMemory(
      BipartiteGraph g,
      int k,
      AlgorithmConfig mhParameters,
      int thresholdForStepThree,
      double thresholdForStepOne,
      boolean retainTopK,
      int topKParameter,
      PrintWriter diagnosticsWriter,
      boolean normalizeColumnsRight,
      boolean normalizeRowsLeft,
      boolean squareWeights,
      String initType,
      boolean applyStepFour,
      double thresholdForStepFour,
      boolean applyLogTransform
  ) {
    //Perform checks on parameters
    if (k != mhParameters.k) {
      throw new RuntimeException(
          String.format(
              "k = %d does not match the k = %d in the AlgorithmConfig object ",
              k, mhParameters.k
          )
      );
    } else if (g.isWeighted()) {
      throw new RuntimeException(
          String.format(
              "Input bipartite graph is weighted -- weighted graphs not supported!"
          )
      );
    }
    this.g = g;
    this.k = k;
    this.mhParameters = mhParameters;
    this.thresholdForStepThree = thresholdForStepThree;
    this.thresholdForStepOne = thresholdForStepOne;
    this.retainTopK = retainTopK;
    if (retainTopK) {
      this.topKParameter = topKParameter;
    } else {
      this.topKParameter = this.g.getNumRightVertices();
    }
    this.diagnosticsWriter = diagnosticsWriter;
    this.squareWeights = squareWeights;
    this.normalizeColumnsRight = normalizeColumnsRight;
    this.normalizeRowsLeft = normalizeRowsLeft;
    this.initType = initType;
    this.applyStepFour = applyStepFour;
    this.thresholdForStepFour = thresholdForStepFour;
    this.applyLogTransform = applyLogTransform;
  }

  /**
   * This method performs Step 1 and Step 2 of SimClusters
   *
   * @return returns a SparseRealMatrix object that contains the latent vector representations for
   * the n vertices on the right, hence an n X k' sparse real matrix. Note that k' might not be
   * equal to k, and the caller must be aware of this.
   */
  private SparseRealMatrix getRepresentationsForRight() {
    //Project the bipartite graph using BipartiteGraph.projectRight() to get
    //an undirected weighted graph on n (number of right vertices) vertices.
    //This is step 1 of SimClusters.
    System.out.println("Computing similarity graph (projection)...");
    Graph projectedGraph = this.g.projectRight(this.thresholdForStepOne);
    //Print statistics of the graph obtained after projection
    System.out.println(
        String.format(
           "Obtained similarity graph (projection) with %d edges and %d vertices",
            projectedGraph.getNumEdges(), projectedGraph.getNumVertices()
        )
    );
    System.out.println(
        String.format("The density of the graph is %f", projectedGraph.getDensity())
    );
    System.out.println(
        "Average edge weight (similarity score) in similarity graph: "
            + projectedGraph.getGlobalAvgWeight()
            + " with standard deviation " + projectedGraph.getWeightStandardDeviation()

    );
    //If squaring weights is specified, then perform it
    if (this.squareWeights) {
      System.out.println("Squaring edge weights in projected graph...");
      projectedGraph.squareWeights();
      System.out.println(
          String.format("The density after squaring is %f", projectedGraph.getDensity())
      );
      System.out.println(
          "Average edge weight of graph after squaring: "
              + projectedGraph.getGlobalAvgWeight()
              + " with standard deviation " + projectedGraph.getWeightStandardDeviation()
      );
    } else {
      System.out.println("No squaring will be performed..");
    }
    //If top-k retention is required then perform it
    if (this.retainTopK) {
      System.out.println("Retaining only top-k edges for every vertex...");
      projectedGraph = projectedGraph.retainTopKEdgesInGraph(this.topKParameter);
      System.out.println(
        String.format(
           "Number of edges after top-k retention: %d edges", projectedGraph.getNumEdges())
      );
      System.out.println(
          String.format("The density after top-k retention is %f", projectedGraph.getDensity())
      );
      System.out.println(
          "Average edge weight after top-k retention: "
              + projectedGraph.getGlobalAvgWeight()
              + " with standard deviation " + projectedGraph.getWeightStandardDeviation()

    );
    } else {
      System.out.println("Retaining all edges in projected graph..");
    }
    int n = projectedGraph.getNumVertices();
    SparseBinaryMatrix z = new SparseBinaryMatrix(n, this.mhParameters.k);
    PrintWriter err = new PrintWriter(System.out);
    // Compute right projection with appropriate threshold
    // If retainTopK is enabled, then retain only top-k edges for every vertex
    if ("randomNeighborhoods".equals(this.initType)) {
      System.out.println("Initializing matrix for step two from random neighborhoods");
      // set allowOverlap = false
      z.initFromBestNeighborhoods(
          projectedGraph, (gp, i) -> this.mhParameters.rng.nextDouble(), false, err
      );
    } else if ("bestNeighborhoods".equals(this.initType)) {
      System.out.println("Initializing matrix from best sub-neighborhoods in terms of conductance");
      z.initFromColSets(
          MHAlgorithm.getRandomBestConductanceSubNeighborhoods(
              projectedGraph, z.getNumCols(), this.mhParameters.rng
          )
      );
    } else if ("nonoverlappingNeighborhoods".equals(this.initType)) {
      System.out.println(
          "Initializing from best non-overlapping sub-neighborhoods in terms of conductance");
      z.initFromColSets(
          MHAlgorithm.getNonOverlappingBestSubNeighborhoods(
              projectedGraph, z.getNumCols(), this.mhParameters.rng
          )
      );
    } else if ("test".equals(this.initType)) {
      System.out.println("Unit test initialization will be used...");
      z.initFromBestNeighborhoods(projectedGraph, Graph::getNeighborhoodConductance,
          false, err);
    } else {
      System.out.println("Initializing factors randomly with 1 nnz/vertex");
      z.initEmptyRowsRandomly(this.mhParameters.rng);
      System.out.println("Random initialization: done");
    }

    //Perform step 2 of the SimClusters algorithm. This involves creating an MHAlgorithm object
    //and calling the optimize() method which runs the Metropolis-Hastings sampling based
    //community detection algorithm on projectedGraph using z as the initialization.
    MHAlgorithm stepTwo =
        new MHAlgorithm(this.mhParameters, projectedGraph, z, this.diagnosticsWriter);
    //Get the cluster assignment
    System.out.println("Running MHAlgorithm on the similarity graph..");
    SparseBinaryMatrix resultOfStepTwo = stepTwo.optimize();
    //Convert and return cluster assignment to a vector containing cluster affiliation scores.
    SparseRealMatrix representationsForRight =
        MHAlgorithm.heuristicallyScoreClusterAssignments(projectedGraph, resultOfStepTwo);
    //Perform column-wise normalization on representationsForRight if specified
    if (this.normalizeColumnsRight) {
      System.out.println("Normalizing columns of right representations..");
      representationsForRight.normalizeToUnitColumn();
    }
    return representationsForRight;
  }

  /**
   * This method computes and returns the latent representation vectors for the m left vertices
   * using the representationsForRight matrix.
   * This is step three of SimClusters algorithm. The thresholding parameter thresholdForStepThree
   * is used for thresholding as mentioned in the class description.
   *
   * @return a SparseRealMatrix object that represents an m X k' matrix containing the latent
   * representations vectors for the m left vertices. Note that k may not be able equal to k'.
   *
   */
  private SparseRealMatrix getRepresentationsForLeft(SparseRealMatrix representationsForRight) {
    //For the rest of the comments in this function k' will be denoted by kp.
    int kp = representationsForRight.getNumCols();
    System.out.println(
        "Expected k: " + this.k + ", Actual k: " + kp
    );
    //We now right-multiply the adjacency matrix A (m X n) of the bipartite graph by the matrix
    //representationsForRight (n X kp) to obtain the representation vectors for the left vertices.
    //We do this without explicitly converting the BipartiteGraph object into a SparseRealMatrix or
    //SparseBinaryMatrix object since the BipartiteGraph object maintains both row-wise and
    //column-wise adjacency information as list of non-zero entries.
    //For the sake of efficiency, instead of computing the dot product of every row of the adjacency
    //matrix with every column of the representationsForRight matrix in order to compute the final
    //product (of A and representationsForRight), we only compute the non-zero entries of the
    //final matrix. This is done by using a hash table (called columnsOfProductMatrix below)
    //which uses (l,i) as key, where l is a column index in the final matrix and i is a row index
    //in the final matrix.
    //To update the hash table, we go over every possible column index j of A, and then use the
    //non-zero entries in the column j of A and the non-zero values in row j of
    //representationsForRight. See below for more details.
    //As mentioned above, we also incorporate a thresholding operation in this step. The idea is to
    //use a Hash table with key (l,i), where l is a column index in the final matrix and i is
    //is a row index in the final matrix, that keeps track of how many right vertices that have a
    //non-zero value in index l of their latent representation (i.e., they "belong" to community l)
    //does the left vertex indexed by i connect to. If this value is less than thresholdForStepThree
    //then we will set the (i,l) entry in the final matrix to zero. The key (l,i) in this Hash table
    //can be thought of as the number of edges vertex i on left has to community l.
    //Hash table for non-zero intersections. This is represented as a nested hash table for the sake
    //of simplicity.
    HashMap<Integer, HashMap<Integer, Double>> columnsOfProductMatrix = new HashMap<>();
    //Hash table for storing number of edges vertex i on left has to community l. To be used for
    //thresholding later.
    HashMap<Integer, HashMap<Integer, Integer>> numEdgesFromLeftToCommunities = new HashMap<>();
    //n = number of columns in A/number of vertices on the right in g
    Integer n = this.g.getNumRightVertices();
    //Go over all column indices
    for (int j = 0; j < n; j++) {
      //Get the left-neighbors for vertex with id j on the right
      //This is basically accessing column j in the matrix A
      int[] neighborsOfVertexJ = this.g.getNeighborsForRightById(j);
      //Get the support of row j in the representationsForRight matrix
      //This is the same as the support of row j in representationsForRight matrix.
      int[] nonZerosColsInRowJ = representationsForRight.getColIdsForRow(j);
      //Get the non-zero values in row j (corresponding to the support)
      double[] nonZeroValuesInRowJ = representationsForRight.getValuesForRow(j);
      //Go over all the elements in the cartesian product: nonZeroColsInRowJ X neighborsOfVertexJ
      //For every (l,i) in this set, update the hash map columnsOfProductMatrix using (l,i) as key:
      //columnsOfProductMatrix[(l,i)] = 0, if (l,i) is not in the map
      //columnsOfProductMatrix[(l,i)]+= representationsForRight[j][l] if (l,i) is in the map
      //update the numEdgesFromLeftToCommunities hash table in a similar manner.
      for (int i : neighborsOfVertexJ) {
        for (int it = 0; it < nonZerosColsInRowJ.length; it++) {
          //compute l as defined
          int l = nonZerosColsInRowJ[it];
          //update the hast table with key (l,i) as mentioned above
          HashMap<Integer, Double> rowIinProductMatrix
              = columnsOfProductMatrix.getOrDefault(l, new HashMap<Integer, Double>());
          double curValOfMatrixAtIL = rowIinProductMatrix.getOrDefault(i, 0.0);
          rowIinProductMatrix.put(i, curValOfMatrixAtIL + nonZeroValuesInRowJ[it]);
          columnsOfProductMatrix.put(l, rowIinProductMatrix);
          //update the number of edges from vertex i on left to community l in a similar manner
          HashMap<Integer, Integer> numEdgesFromVerticesToCommunityL
              = numEdgesFromLeftToCommunities.getOrDefault(l, new HashMap<Integer, Integer>());
          int curNumEdgesFromVertexIToCommunityL
              = numEdgesFromVerticesToCommunityL.getOrDefault(i, 0);
          numEdgesFromVerticesToCommunityL.put(i, curNumEdgesFromVertexIToCommunityL + 1);
          numEdgesFromLeftToCommunities.put(l, numEdgesFromVerticesToCommunityL);
        }
      }
    }
    //We now remove keys (l,i) in the hash table columnsOfProductMatrix for which the number
    //of edges from vertex i (left) to community l is less than the threshold parameter using the
    //information stored in numEdgesFromLeftToCommunities hash table.
    for (int l: numEdgesFromLeftToCommunities.keySet()) {
      HashMap<Integer, Integer> numEdgesFromVerticesToCommunityL
          = numEdgesFromLeftToCommunities.get(l);
      HashMap<Integer, Double> columnLOfProductMatrix = columnsOfProductMatrix.get(l);
      for (int i: numEdgesFromVerticesToCommunityL.keySet()) {
        int curNumEdgesFromVertexIToCommunityL
            = numEdgesFromVerticesToCommunityL.get(i);
        if (curNumEdgesFromVertexIToCommunityL < this.thresholdForStepThree) {
          columnLOfProductMatrix.remove(i);
        }
      }
      columnsOfProductMatrix.put(l, columnLOfProductMatrix);
    }
    //We now use columnsOfProductMatrix to construct a SparseRealMatrix object to represent the
    //m X kp matrix that contains latent representation vectors for the m left vertices. We call it
    //returnMatrix. We need columns as IntSets to define the support of the matrix to be computed
    IntSet[] columnsOfReturnMatrix = new IntSet[kp];
    //Update support of columns of the return matrix using the columnsOfProductMatrix hash map.
    for (int l = 0; l < kp; l++) {
      if (!columnsOfProductMatrix.containsKey(l)) {
        columnsOfReturnMatrix[l] = new IntOpenHashSet(0);
      } else {
        Set<Integer> curCol = columnsOfProductMatrix.get(l).keySet();
        columnsOfReturnMatrix[l] = new IntOpenHashSet(curCol.size());
        for (int rowIndex : curCol) {
          columnsOfReturnMatrix[l].add(rowIndex);
        }
      }
    }
    //m = number of left vertices
    Integer m = this.g.getNumLeftVertices();
    //We now use columnsOfReturnMatrix to define a SparseBinaryMatrix
    SparseBinaryMatrix supportOfReturnMatrix = new SparseBinaryMatrix(m, kp);
    supportOfReturnMatrix.initFromColSets(columnsOfReturnMatrix);
    //Need a 2 dimensional array to store the non-zero double values in the m rows of the return
    //matrix
    double[][] valuesOfReturnMatrix = new double[m][];
    //We populate this 2 dimensional array by traversing the rows of supportOfReturnMatrix.
    for (int row = 0; row < m; row++) {
      int[] nonZeroColsInRow = supportOfReturnMatrix.getRow(row);
      valuesOfReturnMatrix[row] = new double[nonZeroColsInRow.length];
      for (int index = 0; index < nonZeroColsInRow.length; index++) {
        int colId = nonZeroColsInRow[index];
        valuesOfReturnMatrix[row][index] = columnsOfProductMatrix.get(colId).get(row);
      }
    }
    //We can now initialize a SparseRealMatrix using supportOfReturnMatrix and valuesOfReturnMatrix
    SparseRealMatrix returnMatrix =
        new SparseRealMatrix(supportOfReturnMatrix, valuesOfReturnMatrix);
    return returnMatrix;
  }

  /**
   * This method computes and returns the updated latent representation vectors for the n right
   * vertices using the representationsForLeft matrix.
   * The way this is done is by right multiplying the transpose of the m X n bi-adjacency matrix of
   * the bipartite graph with the representationForLeft matrix to obtain an n X k' matrix. The entry
   * in row i and column j of this matrix is then normalized by the value::::
   * sqrt(degree Of Vertex i) X L2Norm(column j of the representationsOfLeft matrix),
   * in other words, the row i and column j entry of the resulting  matrix contains the cosine
   * between the m-dimensional incidence vector of vertex i and column j of the
   * representationsForLeft matrix.
   *
   * This is step four of SimClusters algorithm. The thresholding parameter thresholdForStepFour
   * is also used for thresholding as mentioned in the class description.
   *
   * @return a SparseRealMatrix object that represents an n X k' matrix containing the latent
   * representations vectors for the n left vertices. Note that k may not be able equal to k'.
   *
   * NOTE: We use kp and k' interchangeably in the comments.
   */
  private SparseRealMatrix updateRepresentationsForRight(SparseRealMatrix representationsForLeft) {
    //A nested HashMap to store columns of the nXk' matrix that is obtained by right-multiplying the
    //transpose of the m X n bi-adjacency matrix of the bipartite graph with the m X k' left
    //representations matrix V. Thus, the first key of the HashMap is the column index 0<=j<k', and
    //the keys for the inner HashMap are row-indices 0<=i<n.
    HashMap<Integer, HashMap<Integer, Double>> columnsOfProductMatrix = new HashMap<>();
    //Number of vertices on the left in the bipartite graph g
    Integer m = g.getNumLeftVertices();
    //Number of columns in the product matrix, i.e. k' = the number of columns in the
    //representationsForLeft matrix
    Integer kp = representationsForLeft.getNumCols();
    //Go over all left vertex indices
    for (int l = 0; l < m; l++) {
      //Get the support and non-zero values row l of the leftRepresentations matrix.
      int[] rowLRepLeftSupport = representationsForLeft.getColIdsForRow(l);
      double[] rowLRepLeftValues = representationsForLeft.getValuesForRow(l);
      //Get the (right) neighbors of the vertex L in the bipartite graph
      int[] neighborsOfVertexL = g.getNeighborsForLeftById(l);
      //Go over all the elements in the cartesian product: rowLRepLeftSupport X neighborsOfVertexL
      //For every (j,i) in this set, update the hash map columnsOfProductMatrix using (j,i) as key:
      //columnsOfProductMatrix[(j,i)] = 0, if (j,i) is not in the HashMap as a key yet,
      //columnsOfProductMatrix[(j,i)]+= representationsForLeft[l][j] if (j,i) is present as a key.
      for (int it = 0; it < rowLRepLeftSupport.length; it++) {
        //compute j as defined
        int j = rowLRepLeftSupport[it];
        //Get column j of the product matrix from the HashMap
        HashMap<Integer, Double> colJInProductMatrix
              = columnsOfProductMatrix.getOrDefault(j, new HashMap<Integer, Double>());
        for (int i : neighborsOfVertexL) {
          //update the hash table with key (j,i) as mentioned above
          double curValOfMatrixAtColJRowI = colJInProductMatrix.getOrDefault(i, 0.0);
          colJInProductMatrix.put(i, curValOfMatrixAtColJRowI + rowLRepLeftValues[it]);
        }
        //Update column j of the product matrix in the HashMap
        columnsOfProductMatrix.put(j, colJInProductMatrix);
      }
    }
    //We now update the value for each key (j,i) in the hash table columnsOfProductMatrix to:
    //columnOfProductMatrix[(j,i)]/Sqrt(NumNeighborsOfVertexI) X L2Norm(ColJofRepresentationForLeft
    //Additionally, if this value is less than thresholdForStepFour then we remove it from the table
    //In case of updating the existing columnsOfProductMatrix HashMap, we create a new HashMap.
    HashMap<Integer, HashMap<Integer, Double>> columnsOfProductMatrixUpdated =
        new HashMap<>();
    for (int j: columnsOfProductMatrix.keySet()) {
      HashMap<Integer, Double> columnJOfProductMatrix = columnsOfProductMatrix.get(j);
      HashMap<Integer, Double> newColumnJOfProductMatrix = new HashMap<>();
      double l2NormOfColumnJinLeftRepMatrix = representationsForLeft.getColumnNorm(j);
      for (int i: columnJOfProductMatrix.keySet()) {
        double sqrtOfDegreeOfVertexI = Math.sqrt(g.getNeighborsForRightById(i).length);
        double normalizationFactor =
            l2NormOfColumnJinLeftRepMatrix * sqrtOfDegreeOfVertexI;
        double newValueOfProductMatrixAtColJRowI =
            columnJOfProductMatrix.get(i) / normalizationFactor;
        //Remove value if it is less than the
        if (newValueOfProductMatrixAtColJRowI - thresholdForStepFour > 1e-5) {
          newColumnJOfProductMatrix.put(i, newValueOfProductMatrixAtColJRowI);
        }
      }
      //if newColumnJofProductMatrix has no keys, i.e. it's empty, then don't add it to new table
      if (newColumnJOfProductMatrix.size() > 0) {
        columnsOfProductMatrixUpdated.put(j, newColumnJOfProductMatrix);
      }
    }
    //Set columnsOfProductMatrix to columnsOfProductMatrixUpdated
    columnsOfProductMatrix = columnsOfProductMatrixUpdated;
    //We now use columnsOfProductMatrix to construct a SparseRealMatrix object to represent the
    //n X k' matrix that contains latent representation vectors for the n right vertices. We call it
    //returnMatrix. We need columns as IntSets to define the support of the matrix to be computed
    IntSet[] columnsOfReturnMatrix = new IntSet[kp];
    //Update support of columns of the return matrix using the columnsOfProductMatrix hash map.
    for (int j = 0; j < kp; j++) {
      if (!columnsOfProductMatrix.containsKey(j)) {
        columnsOfReturnMatrix[j] = new IntOpenHashSet(0);
      } else {
        Set<Integer> curCol = columnsOfProductMatrix.get(j).keySet();
        columnsOfReturnMatrix[j] = new IntOpenHashSet(curCol.size());
        for (int rowIndex : curCol) {
          columnsOfReturnMatrix[j].add(rowIndex);
        }
      }
    }
    //n = number of right vertices
    Integer n = this.g.getNumRightVertices();
    //We now use columnsOfReturnMatrix to define a SparseBinaryMatrix
    SparseBinaryMatrix supportOfReturnMatrix = new SparseBinaryMatrix(n, kp);
    supportOfReturnMatrix.initFromColSets(columnsOfReturnMatrix);
    //Need a 2 dimensional array to store the non-zero double values in the n rows of the return
    //matrix
    double[][] valuesOfReturnMatrix = new double[n][];
    //We populate this 2 dimensional array by traversing the rows of supportOfReturnMatrix.
    for (int row = 0; row < n; row++) {
      int[] nonZeroColsInRow = supportOfReturnMatrix.getRow(row);
      valuesOfReturnMatrix[row] = new double[nonZeroColsInRow.length];
      for (int index = 0; index < nonZeroColsInRow.length; index++) {
        int colId = nonZeroColsInRow[index];
        valuesOfReturnMatrix[row][index] = columnsOfProductMatrix.get(colId).get(row);
      }
    }
    //We can now initialize a SparseRealMatrix using supportOfReturnMatrix and valuesOfReturnMatrix
    SparseRealMatrix returnMatrix =
        new SparseRealMatrix(supportOfReturnMatrix, valuesOfReturnMatrix);
    return returnMatrix;
  }

  /**
   * Calling this method performs steps one, two, and three of the SimClusters in-memory algorithm
   * in succession on the input bipartite graph and returns the latent representation vectors
   * (of dimension k) for all the vertices.
   * Additionally, step 4 is performed if the applyStepFour variable is set to true.
   *
   * @return returns the representations of the left and right vertices as matrices bundled together
   * in a BipartiteSparseRepresentations object.
   */
  public BipartiteSparseRepresentations runSimClusters() {
    //Perform steps one and two to get right representations
    SparseRealMatrix representationsForRight = this.getRepresentationsForRight();

    //Perform step three to get the representation vectors for the left vertices using the
    //representationsForRight matrix. This will be an m X k' matrix.
    SparseRealMatrix representationsForLeft =
        this.getRepresentationsForLeft(representationsForRight);
    //If rows have to normalized in the left matrix, do that
    if (this.normalizeRowsLeft) {
      System.out.println("Normalizing left embeddings");
      representationsForLeft.normalizeToUnitColumn();
    }
    //If applyLogTransform is true, then apply the log10(1+x) transform
    if (this.applyLogTransform) {
      System.out.println("Applying coordinate-wise log10(1+x) transform");
      representationsForLeft.applyCoordinateWiseLogTransform();
    }
    //If fourth step needs to performed
    if (this.applyStepFour) {
      System.out.println("Applying Step Four");
      double sparsityBeforeStepFour = representationsForRight.getAverageNNZ();
      System.out.println(
        String.format(
            "Average NNZ per row in right representation matrix before step 4: %f",
            sparsityBeforeStepFour
        )
      );
      representationsForRight = this.updateRepresentationsForRight(representationsForLeft);
      double sparsityAfterStepFour = representationsForRight.getAverageNNZ();
      System.out.println(
        String.format(
            "Average NNZ per row in right representation matrix after step 4: %f",
            sparsityAfterStepFour
        )
      );
    }
    //We now return a BipartiteSparseRepresentation object that contains the representations for the
    //left and right vertices.
    BipartiteSparseRepresentations newRepresentation =
        new BipartiteSparseRepresentations(
            representationsForLeft,
            representationsForRight
        );
    return newRepresentation;
  }

}
