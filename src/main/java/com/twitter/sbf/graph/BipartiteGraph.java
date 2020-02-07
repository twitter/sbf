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

package com.twitter.sbf.graph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;

import com.twitter.sbf.util.FileLineIterator;
import com.twitter.sbf.util.SimpleIterator;


/**
 * Class for representing unweighted bipartite graph as adjacency list (sets)
 * The vertices on the left hand side are 0-indexed
 * The vertices on the right hand side are also 0-indexed, separate from the ones on the left
 */
public class BipartiteGraph {

  //Number of edges, left vertices, and right vertices
  private long numEdges;
  private int numLeftVertices;
  private int numRightVertices;
  //Adjacency lists for each vertex, neighborsForLeft are the adjacency lists for the left vertices
  //similarly for right vertices
  //Each adjacency list is sorted for efficient dot product
  private int[][] neighborsForLeft;
  private int[][] neighborsForRight;

  //Constructor
  public BipartiteGraph(
      long numEdges, int numLeftVertices, int numRightVertices,
      int[][] neighborsForLeft, int[][] neighborsForRight
  ) {
    this.numEdges = numEdges;
    this.numLeftVertices = numLeftVertices;
    this.numRightVertices = numRightVertices;
    this.neighborsForLeft = neighborsForLeft;
    this.neighborsForRight = neighborsForRight;
    this.areListsSortedAndSymmetric();
  }

  //access methods
  public long getNumEdges() {
    return this.numEdges;
  }

  public int getNumLeftVertices() {
    return this.numLeftVertices;
  }

  public int getNumRightVertices() {
    return this.numRightVertices;
  }

  /**
   * Returns the adjacency list for vertex 'id' on the left side
   */
  public int[] getNeighborsForLeftById(int id) {
    if (id < this.numLeftVertices) {
      return this.neighborsForLeft[id];
    } else {
      throw new RuntimeException(String.format("No left side vertex with id %d", id));
    }
  }

  /**
   * Returns the adjacency list for vertex 'id' on the right side
   */
  public int[] getNeighborsForRightById(int id) {
    if (id < this.numRightVertices) {
      return this.neighborsForRight[id];
    } else {
      throw new RuntimeException(String.format("No right side vertex with id %d", id));
    }
  }

  //get edge density of graph
  public double getEdgeDensity() {
    return ((double) this.numEdges) / ((long) this.numRightVertices * this.numLeftVertices);
  }

  /**
   * Constructor
   * The input must be in a modified METIS format for bipartite graphs:
   * First line contains numLeftVertices numRightVertices numEdges
   * The next numLeftVertices lines contain adjacency lists for left side vertices (1-indexed)
   * The next numrightVertices lines contain adjacency lists for right side vertices (1-indexed)
   * There could be comment lines in between that begin with a '%' which will ignored.
   */
  public BipartiteGraph(SimpleIterator<String> lines) {
    Optional<String> lineOpt;

    lineOpt = lines.next();
    if (!lineOpt.isPresent()) {
      throw new RuntimeException("No lines in input!");
    } else {
      String[] tokens = lineOpt.get().trim().split("\\s+");
      assert tokens.length == 3;
      this.numLeftVertices = Integer.parseInt(tokens[0]);
      this.numRightVertices = Integer.parseInt(tokens[1]);
      //this is the number of vertices we expect there to be in the graph
      //we will ensure that the number of edges in the string matches in this.
      long expectedNumEdges = Long.parseLong(tokens[2]);
      //initialize the adjacency lists for vertices on the left and the right
      this.neighborsForLeft = new int[numLeftVertices][];
      this.neighborsForRight = new int[numRightVertices][];

      //temp variables and loop to read in graph data from the string
      //checks incorporated to make sure input is valid and consistent
      int nodeId = 0;
      long edgesSoFar = 0;
      while (true) {
        lineOpt = lines.next();
        if (!lineOpt.isPresent()) {
          break;
        } else {
          int[] myNeighbors;
          String line = lineOpt.get();
          if (line.startsWith("%")) {
            continue;
          }
          line = line.trim();
          // Fields are separated by spaces
          tokens = line.split("\\s+");
          if (line.isEmpty()) { // vertices with no neighbor
            myNeighbors = new int[0];
          } else { // The rest are adjacency lists
            int numNeighbors = tokens.length;
            edgesSoFar += numNeighbors;
            myNeighbors = new int[numNeighbors];
            for (int i = 0; i < numNeighbors; i++) {
              int nId = Integer.parseInt(tokens[i]) - 1; // 0-based indexing
              if (nId < 0) {
                throw new IllegalArgumentException(
                    String.format(
                        "Line %d: neighbor %d is less than 1 (all vertex ids are positive)\n",
                        nodeId + 1, nId + 1
                    )
                );
              } else if (nodeId < this.numLeftVertices && nId >= this.numRightVertices) {
                throw new IllegalStateException(
                    String.format(
                        "Line %d: neighbor %d is more than the number of right side vertices %d\n",
                        nodeId + 1, nId + 1, this.numRightVertices
                    )
                );
              } else if (nodeId > this.numLeftVertices && nId >= this.numLeftVertices) {
                throw new IllegalArgumentException(
                    String.format(
                        "Line %d: neighbor %d is more than the number of left side vertices %d\n",
                        nodeId + 1, nId + 1, this.numLeftVertices
                    )
                );
              }
              myNeighbors[i] = nId;
            }
            Arrays.sort(myNeighbors);
          }
          if (nodeId < this.numLeftVertices) {
            this.neighborsForLeft[nodeId++] = myNeighbors;
          } else if (nodeId < this.numLeftVertices + this.numRightVertices) {
            this.neighborsForRight[(nodeId++) - this.numLeftVertices] = myNeighbors;
          } else {
            int totalVertices = this.numLeftVertices + this.numRightVertices;
            throw new IllegalArgumentException("Expected number of nodes " + totalVertices
                + " is inconsistent with number of adjacency lists");
          }
          if (nodeId % 100000 == 0) {
            System.err.print("\rDone loading " + nodeId + " lines.");
          }
        }
      }
      if (edgesSoFar % 2 != 0) {
        throw new RuntimeException("Got odd number of edges: "
            + edgesSoFar + ". Graph is not valid undirected graph.");
      } else if (edgesSoFar != expectedNumEdges * 2) {
        throw new RuntimeException("Expected " + 2 * expectedNumEdges
            + " (counting both directions), but got " + edgesSoFar + "  edges.");
      } else if (nodeId < this.numLeftVertices + this.numRightVertices) {
        throw new IllegalArgumentException("Incomplete METIS data file!");
      } else {
        this.numEdges = expectedNumEdges;
      }
    }
    this.areListsSortedAndSymmetric();
  }

  /**
   * Constructor
   * Input file must contain the bipartite graph in METIS format
   * (See comments above previous function)
   *
   * @param filename (required)
   */
  public static BipartiteGraph fromFile(String filename) throws IOException {
    BipartiteGraph ret = null;
    BufferedReader br = new BufferedReader(new FileReader(filename));
    SimpleIterator<String> lines = new FileLineIterator(br);
    ret = new BipartiteGraph(lines);
    br.close();
    return ret;
  }

  /**
   * Method to return the transpose of the bipartite graph, i.e. swap left and right side
   *
   * @return a new BipartiteGraph object containing the transposed version of the calling
   * Bipartite graph.
   */
  public BipartiteGraph getTranspose() {
    return new BipartiteGraph(
        this.numEdges, this.numRightVertices,
        this.numLeftVertices, this.neighborsForRight, this.neighborsForLeft
    );
  }

  /**
   * Method to project the bipartite graph onto the right side using a weighted projection.
   * The result is a weighted undirected graph on n = numRightVertices vertices.
   * The weight of the edge between i and j (i and j are vertices on the right side) is
   * |N(i) intersection N(j)|/Sqrt(|N(i)|*|N(j)|)
   * where N(i) and N(j) are the sets of neighbors of i and j on the left hand side respectively.
   * Also supports thresholding: if the weight of an edge in the projected graph is less than
   * the user specified threshold, the edge is deleted. Tolerance parameter is set to 1e-5.
   *
   * @param threshold
   */
  public Graph projectRight(double threshold) {
    /*
     * Maintain a hash table T where the keys are (i,j) for 1 <= i < j <= numRightVertices,
     * and the values are integer counts whose final value will represent
     * the size of the intersection of the adjacency lists of vertices i and j on the right side,
     * i.e. |N(i) intersection N(j)|.
     * To populate the hash table go through every vertex on the left side, and
     * for every such vertex, look at all the pairs of vertices (i,j), i<j, formed
     * using vertices i and j in its adjacency list.
     * For every pair (i,j), update the hash table as follows: T[(i,j)] = T[(i,j)] + 1, if
     * (i,j) is already present in the table, otherwise set T[(i,j)] = 1.
     * At the end of the process, if N(i) and N(j) don't intersect, then (i,j) is not in T,
     * otherwise T[(i,j)] is exactly |N(i) intersection N(j)|.
     * We can now compute the adjacency lists and edge-weight lists for every vertex
     * in the projected graph using the table T[(i,j)]. We also perform thresholding at this stage
     * and delete edges whose weights are below the threshold.
     * Since the Graph class constructor requires that we supply adjacency and weight lists
     * in sorted order, we incur some additional cost in generating sorting lists. This turns out
     * to be an additional multiplicative factor equaling log(numRightVertices)<<log(10^8)~30.
     * Also, we implement the hash table T as a nested two-level HashTable for efficiency reasons.
     */
    //Define and populate the hash table T (we refer to it as pairsOnTheRightThatShare in the code)
    //Iterate through the vertices on the left side
    HashMap<Integer, HashMap<Integer, Integer>> pairsOnTheRightThatShare = new HashMap<>();
    for (int leftVertex = 0; leftVertex < this.numLeftVertices; leftVertex++) {
      //Get the adjacency list of the current left vertex
      int[] curAdjList = this.neighborsForLeft[leftVertex];
      //Iterate over all pairs i<j in the adjacency list of the current left vertex we are on
      for (int it1 = 0; it1 < curAdjList.length; it1++) {
        for (int it2 = it1 + 1; it2 < curAdjList.length; it2++) {
          //Define i and j as above.
          //Note that since the adjacency lists in the bipartite graph are sorted,
          //we have that i < j.
          int i = curAdjList[it1];
          int j = curAdjList[it2];
          HashMap<Integer, Integer> mapOfI
              = pairsOnTheRightThatShare.getOrDefault(i, new HashMap<>());
          Integer curValOfIAndJ = mapOfI.getOrDefault(j, 0);
          mapOfI.put(j, curValOfIAndJ + 1);
          pairsOnTheRightThatShare.put(i, mapOfI);
          //For convenience sake, we also maintain the same count at T[(j,i)].
          //This symmetry allows us to construct the adjacency lists of the projected graph
          //with ease in the next step.
          HashMap<Integer, Integer> mapOfJ =
              pairsOnTheRightThatShare.getOrDefault(j, new HashMap<>());
          mapOfJ.put(i, curValOfIAndJ + 1);
          pairsOnTheRightThatShare.put(j, mapOfJ);
        }
      }
    }
    //Use the hash table T to create adjacency lists and weight lists for every vertex
    //in the projected graph. Use the user specified threshold along with a tolerance of 1e-5 to
    //delete edges with weight below the threshold.
    //We also maintain the number of edges in the projected graph.
    long numEdgesProjected = 0L;
    int[][] adjListsProjected = new int[this.numRightVertices][];
    float[][] wtListsProjected = new float[this.numRightVertices][];
    //We also maintain a HashTable to keep track of edge weights that have been
    //already computed to avoid recomputing twice for the same edge.
    HashMap<Integer, HashMap<Integer, Float>> computedEdgeWeights = new HashMap<>();
    for (int i = 0; i < this.numRightVertices; i++) {
      HashMap<Integer, Integer> mapOfI = pairsOnTheRightThatShare.getOrDefault(i, new HashMap<>());
      //Define temporary ArrayLists to store the adjacency list and weights for vertex i
      //in the projected graph after thresholding.
      ArrayList<Integer> thresholdedAdjacencyListForI = new ArrayList<Integer>(0);
      ArrayList<Float> thresholdedWtListForI = new ArrayList<Float>(0);
      //The number of neighbors of vertex i in the projected graph after thresholding
      int thresholdedNumOfNeighborsOfI = 0;
      //If the list of neighbors of i is empty, then there is nothing more to do.
      if (mapOfI.size() > 0) {
        //Convert the set of keys of mapOfI into an array. These are neighbors of i
        //in the projected graph pre-thresholding.
        Integer[] keysOfI = mapOfI.keySet().toArray(new Integer[mapOfI.size()]);
        //Sort the array of neighbors of i by their index/id.
        Arrays.sort(keysOfI);
        //Iterate through the neighbors in sorted list and set the weight and adjacency for each
        for (int nId = 0; nId < mapOfI.size(); nId++) {
          //Get the label of the neighbor
          int j = keysOfI[nId];
          //float variable to hold the weight of the edge between i and j
          float weightOfIJEdge;
          //If j < i then the weight for the edge between i and j in the projected graph
          // would have been computed during the iteration through
          // the neighbors of j and so we can just look it up
          //from the hash table.
          if (j < i) {
            weightOfIJEdge = computedEdgeWeights.get(j).get(i);
          } else {
            //Recall that T[(i,j)] (and T[(j,i)]) now contains |N(i) intersection N(j)|
            int intersectionSize = mapOfI.get(j);
            //To compute weight of the edge between i and j, we must normalize
            //the size of the intersection by Sqrt(degreeOfI * degreeOfJ),
            //where degreeOfI and degreeOfJ are the degrees of i and j in the bipartite graph.
            int degreeOfI = this.getNeighborsForRightById(i).length;
            int degreeOfJ = this.getNeighborsForRightById(j).length;
            weightOfIJEdge = intersectionSize / (float) Math.sqrt(degreeOfI * degreeOfJ);
            //Add the computed weight to the hash table
            HashMap<Integer, Float> tempWeightMap
                = computedEdgeWeights.getOrDefault(i, new HashMap<Integer, Float>());
            tempWeightMap.put(j, weightOfIJEdge);
            computedEdgeWeights.put(i, tempWeightMap);
          }
          //Update the adjacency and weight list of vertex i IF the weight of the edge
          //between i and j is greater than the threshold (using tolerance 1e-5)
          if (weightOfIJEdge - threshold > 1e-5) {
            thresholdedNumOfNeighborsOfI++;
            thresholdedAdjacencyListForI.add(j);
            thresholdedWtListForI.add(weightOfIJEdge);
          }
        }
      }
      //Update number of edges. Note this is double counting. We will divide by two later.
      numEdgesProjected += thresholdedNumOfNeighborsOfI;
      //Initialize the adjacency and weight list of vertex i in the projected graph
      adjListsProjected[i] = new int[thresholdedNumOfNeighborsOfI];
      wtListsProjected[i] = new float[thresholdedNumOfNeighborsOfI];
      //Convert the ArrayLists to primitive arrays
      for (int index = 0; index < thresholdedNumOfNeighborsOfI; index++) {
        adjListsProjected[i][index] = thresholdedAdjacencyListForI.get(index);
        wtListsProjected[i][index] = thresholdedWtListForI.get(index);
      }
    }
    //Every edge has been counted twice, so we divide by two
    //Instantiate a new instance of Graph and return
    return new Graph(
        this.numRightVertices, numEdgesProjected / 2, adjListsProjected, wtListsProjected
    );
  }

  //Method to project the bipartite graph on to the left side (similar to projectRight)
  public Graph projectLeft(double threshold) {
    return this.getTranspose().projectRight(threshold);
  }

  /**
   * Method to check if an array of integers is in sorted order
   * @param list
   * @return value
   */
  private void checkIsSorted(int[] list) throws IllegalArgumentException {
    for (int j = 0; j < list.length - 1; j++) {
        if (list[j] >= list[j + 1]) {
          throw new IllegalArgumentException("Adjacency lists are not strictly increasing! ");
        }
    }
  }

  //Method to check if the adjacency lists are strictly increasing and symmetric
  private void areListsSortedAndSymmetric() throws IllegalArgumentException {
    //Make sure that the adjacency lists for each vertex are in strictly increasing order
    for (int i = 0; i < this.numLeftVertices; i++) {
      int[] curLeftList = this.neighborsForLeft[i];
      checkIsSorted(curLeftList);
    }
    for (int i = 0; i < this.numRightVertices; i++) {
      int[] curRightList = this.neighborsForRight[i];
      checkIsSorted(curRightList);
    }

    //Make sure that the bipartite graph is symmetric,
    // i.e. if i is the list of j, then j is also in the list of i
    //Use HashSets for doing the checks
    ArrayList<HashSet<Integer>> neighborsForLeftHashed =
        new ArrayList<HashSet<Integer>>(this.numLeftVertices);
    ArrayList<HashSet<Integer>> neighborsForRightHashed =
        new ArrayList<HashSet<Integer>>(this.numRightVertices);
    for (int[] curList : this.neighborsForLeft) {
      HashSet<Integer> tempSet = new HashSet<Integer>(0);
      for (int member : curList) {
        tempSet.add(member);
      }
      neighborsForLeftHashed.add(tempSet);
    }
    for (int[] curList : this.neighborsForRight) {
      HashSet<Integer> tempSet = new HashSet<Integer>(0);
      for (int member : curList) {
        tempSet.add(member);
      }
      neighborsForRightHashed.add(tempSet);
    }
    for (int i = 0; i < this.numLeftVertices; i++) {
      int[] curList = this.neighborsForLeft[i];
      for (int member : curList) {
        if (!neighborsForRightHashed.get(member).contains(i)) {
          throw new IllegalArgumentException("Adjacency lists are not symmetric! ");
        }
      }
    }
    for (int i = 0; i < this.numRightVertices; i++) {
      int[] curList = this.neighborsForRight[i];
      for (int member : curList) {
        if (!neighborsForLeftHashed.get(member).contains(i)) {
          throw new IllegalArgumentException("Adjacency lists are not symmetric! ");
        }
      }
    }
  }
  //This method checks if the graph is weighted. Currently we don't have support for weighted
  //bipartite graphs. If in the future the class is extended to weighted bipartite graphs then
  //the method must be modified accordingly:
  // replace "return false" with "return (weights != null)" where weights is a double[][] array
  //that represents the weights of edges emanating from the left vertices.
  public Boolean isWeighted() {
    return false;
  }
}
