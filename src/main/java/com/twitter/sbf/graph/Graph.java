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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.ImmutableMap;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.random.RandomAdaptor;

import com.twitter.sbf.util.FileLineIterator;
import com.twitter.sbf.util.SimpleIterator;

import it.unimi.dsi.fastutil.ints.IntSet;

/**
 * Undirected graph, can be either weighted or not, represented as adjacency lists.
 * Neighbors of a vertex is stored in an ascending array.
 * Vertex ids are 0-based.
 */
public class Graph {
  private int numVertices; // number of vertices
  private long numEdges; // number of edges
  private int[][] neighbors; // neighbors[i] is the ascendingly sorted array of neighbors of i
  private float[][] weights; // same cardinality as neighbors, gives weight of corresponding edge
  private float[] weightedOutDegrees; // same cardinality as neighbors and weights,
  private double globalAvgWeight;
  // ith entry gives sum of outgoing edges for node i.
  private int[] allVertexIds;
  private EnumeratedIntegerDistribution degreeDist = null;

  public Graph(int numVertices, long numEdges, int[][] nbrs, float[][] wts) {
    this.numVertices = numVertices;
    this.numEdges = numEdges;
    this.neighbors = nbrs;
    this.weights = wts;
    if (isWeighted()) {
      double totalSum = 0;
      this.weightedOutDegrees = new float[this.numVertices];
      for (int i = 0; i < this.numVertices; i++) {
        this.weightedOutDegrees[i] = 0;
        for (int j = 0; j < this.weights[i].length; j++) {
          if (Float.isNaN(this.weights[i][j])) {
            throw new RuntimeException(
                "Weight from node " + i + " to node " + j + " (zero-indexed) is NaN!"
            );
          }
          this.weightedOutDegrees[i] += this.weights[i][j];
          totalSum += this.weights[i][j];
        }
      }
      // divide by 2 since each edge is present twice
      this.globalAvgWeight = totalSum / this.numEdges / 2;
    }
    checkSizesOfNbrsAndWts();
    checkSymmetryAndSorting();
  }

  /**
   * Calculate conductance for vertex i in Graph this.
   * @param i
   * @return
   */
  public double getNeighborhoodConductance(int i) {
    if (this.getDegree(i) < 2) {
      return 1.0;
    }
    double cut = 0;
    double vol = 0;

    for (Integer vertexId : this.getNeighbors(i)) {
      Iterator<WeightedId> wni = this.getWeightedNeighborsIterator(vertexId);
      while (wni.hasNext()) {
        WeightedId wn = wni.next();
        if (!this.isNeighbors(wn.neighborId, i) && wn.neighborId != i) {
          cut += wn.weight;
        }
      }
      vol += this.getWeightedOutDegree(vertexId);
    }
    vol += this.getWeightedOutDegree(i);

    double totalEdgeWt = 0;
    if (this.isWeighted()) {
      totalEdgeWt = 2.0 * this.getNumEdges() * this.getGlobalAvgWeight();
    } else {
      totalEdgeWt = 2.0 * this.getNumEdges();
    }

    if (vol != totalEdgeWt && vol != 0) {
      return cut / Math.min(vol, totalEdgeWt - vol);
    } else {
      return 0.0;
    }
  }

  /**
   * Compute conductance of S = (vertex i, neighbors)
   *
   * @param i (required) vertex id
   * @return conductance score of S
   */
  public double getNeighborhoodConductanceUnweighted(int i) {
    if (this.getDegree(i) < 2) {
      return 1.0;
    }
    double cut = 0;
    double vol = 0;
    for (Integer vertexId : this.getNeighbors(i)) {
      for (Integer neighborId : this.getNeighbors(vertexId)) {
        if (!this.isNeighbors(neighborId, i) && neighborId != i) {
          cut += 1;
        }
      }
      vol += this.getDegree(vertexId);
    }
    vol += this.getDegree(i);
    double edge2 = 2.0 * this.getNumEdges();
    if (vol != edge2 && vol != 0) {
      return cut / Math.min(vol, edge2 - vol);
    } else {
      return 0.0;
    }
  }

  public int[][] getNeighbors() {
    return neighbors;
  }

  public boolean isWeighted() {
    return weights != null;
  }

  public double getDensity() {
    return (2.0 * this.numEdges / this.numVertices) / (this.numVertices - 1);
  }

  /**
   * Get average weight of an edge in the graph as a whole.
   * @return average weight of an edge in the graph as a whole. In the case of an unweighted graph,
   * returns 1.
   */
  public double getGlobalAvgWeight() {
    if (isWeighted()) {
      return globalAvgWeight;
    } else {
      return 1;
    }
  }

  /**
   * Get the standard deviation of the edge weights.
   * @return returns the standard deviation if the graph is weighted otherwise returns 0
   */
  public double getWeightStandardDeviation() {
    if (!isWeighted()) {
      return 0;
    }
    //Compute average
    double avgWeight = this.getGlobalAvgWeight();
    //iterate through all edges and sum up the squares of the deviations
    double deviationSquareSums = 0.0;
    for (float[] weightList: this.weights) {
      for (float wt: weightList) {
        deviationSquareSums += Math.pow(wt - avgWeight, 2.0);
      }
    }
    double variance = 0.0;
    if (this.getNumEdges() > 1) {
      variance = deviationSquareSums / (this.getNumEdges() - 1);
    }
    return Math.sqrt(variance);
  }

  public double getGlobalWeight() {
    return getGlobalAvgWeight() * getNumEdges() * 2.0;
  }


  /**
   * Checks if (a) the adjacency matrix of the graph is symmetric, and
   * (b) the neighbors of each node are in ascending order (if the sorting is not maintained,
   * many methods such as getWeightOfEdge() will give incorrect results)
   */
  public void checkSymmetryAndSorting() {
    for (int i = 0; i < this.getNumVertices(); i++) {
      int[] neighborsOfI = this.getNeighbors(i);
      for (int j = 0; j < neighborsOfI.length; j++) {
        if (j > 0 && neighborsOfI[j] <= neighborsOfI[j - 1]) {
          throw new IllegalStateException(
            "Input graph does not have ascendingly sorted neighbors! "
            + String.format(
              "Neighbor no. %d (value %d) of node %d is less than neighbor no. %d (value %d)",
              // add 1s to everything to get back 1-indexing
              j + 1, neighborsOfI[j] + 1, i + 1, j, neighborsOfI[j - 1] + 1
            )
          );
        }
        if (neighborsOfI[j] < i) {
          double diff = Math.abs(
            this.getWeightOfJInI(neighborsOfI[j], i) - this.getWeightOfJInI(i, neighborsOfI[j])
          );
          if (diff > 1e-5) {
            throw new IllegalStateException("Input graph is not symmetric! "
              + String.format(
                "Weight of (%d, %d) is %.2g, while weight of (%d, %d) is %.2g",
                i + 1, neighborsOfI[j] + 1, this.getWeightOfJInI(i, neighborsOfI[j]),
                neighborsOfI[j] + 1, i + 1, this.getWeightOfJInI(neighborsOfI[j], i)
              )
            );
          }
        }
      }
    }
  }

  private void checkSizesOfNbrsAndWts() {
    if (this.neighbors != null) {
      int numNeighbors = 0;
      for (int[] arr : this.neighbors) {
        numNeighbors += arr.length;
      }
      assert numNeighbors == this.numEdges * 2
          : String.format(
              "Expected size of neighbors %d, actual size of neighbors %d",
                this.numEdges * 2, numNeighbors
            );
    }
    if (this.weights != null) {
      int numNeighbors = 0;
      for (float[] arr : this.weights) {
        numNeighbors += arr.length;
      }
      assert numNeighbors == this.numEdges * 2
          : String.format(
              "Expected size of weights %d, actual size of weights %d",
                this.numEdges * 2, numNeighbors
            );
    }
  }

  public Graph(int numVertices, long numEdges, int[][] nbrs) {
    this(numVertices, numEdges, nbrs, null);
  }

  /**
   * Load unweighted graph into nbrs argument
   *
   * @param lines iterator over lines
   * @param nbrs 2d array with neighbors of each vertex. this needs to pre-allocated in the
   * outer dimension and will be filled in by this function.
   * @param expectedNumEdges number of edges expected by caller
   * @return actual number of edges
   */
  private static long loadUnweightedGraph(
      SimpleIterator<String> lines,
      int[][] nbrs,
      long expectedNumEdges) {
    Optional<String> lineOpt;
    int nodeId = 0;
    long edgesSoFar = 0;
    assert nbrs != null;
    while (true) {
      lineOpt = lines.next();
      if (!lineOpt.isPresent()) {
        break;
      } else {
        String line = lineOpt.get();
        if (line.startsWith("%")) {
          continue;
        }
        line = line.trim();
        // Fields are separated by spaces
        String[] tokens = line.split("\\s+");
        if (line.isEmpty()) { // vertices with no neighbor
          nbrs[nodeId++] = new int[0];
        } else { // The rest are adjacency lists
          int numNeighbors = tokens.length;
          edgesSoFar += numNeighbors;
          int[] myNeighbors = new int[numNeighbors];
          for (int i = 0; i < numNeighbors; i++) {
            int nId = Integer.parseInt(tokens[i]) - 1; // 0-based indexing
            if (nId < 0) {
              throw new IllegalStateException(
                  String.format(
                      "Line %d: neighbor %d is less than 1 (all vertex ids are positive)\n",
                      nodeId + 1, nId + 1
                  )
              );
            }
            myNeighbors[i] = nId;
          }
          Arrays.sort(myNeighbors);
          nbrs[nodeId++] = myNeighbors;
        }
        if (nodeId > nbrs.length) {
          throw new RuntimeException("Expected number of nodes " + nbrs.length
              + " is inconsistent with number of adjacency lists");
        }
        if (nodeId % 100000 == 0) {
          System.err.print("\rDone loading " + nodeId + " lines.");
        }
      }
    }
    if (edgesSoFar % 2 != 0) {
      throw new RuntimeException("Got odd number of edges (counting edges from both directions) - "
          + edgesSoFar + ". Graph is not valid undirected graph.");
    }
    if (edgesSoFar != expectedNumEdges * 2) {
      System.err.println("Expected " + 2 * expectedNumEdges
          + " (counting both directions), but got " + edgesSoFar + " many edges.");
    }
    return edgesSoFar / 2;
  }

  static int getIdFromIdWeightString(String s) {
    int separatorIndex = s.indexOf(" ");
    assert separatorIndex >= 0;
    return Integer.parseInt(s.substring(0, separatorIndex));
  }

  private static float getWeightFromIdWeightString(String s) {
    int separatorIndex = s.indexOf(" ");
    assert separatorIndex >= 0;
    return Float.parseFloat(s.substring(separatorIndex + 1, s.length()));
  }

  static class IdWeightStringComparator implements Comparator<String> {
    @Override
    public int compare(String o1, String o2) {
      return Integer.compare(getIdFromIdWeightString(o1), getIdFromIdWeightString(o2));
    }
  }

  /**
   * Get a distribution object representing the degree distribution of the vertices of this graph,
   * suitable for sampling from. Note this uses the unweighted degree distribution.
   * @param rng Random number generator
   * @return distribution object representing the degree distribution of the vertices of this graph,
   * suitable for sampling from.
   */
  public EnumeratedIntegerDistribution getDegreeDistribution(RandomAdaptor rng) {
    if (degreeDist == null) {
      double[] degreesNormalized = new double[numVertices];
      int[] vertices = new int[numVertices];
      for (int i = 0; i < numVertices; i++) {
        degreesNormalized[i] = this.getDegree(i) / (2.0 * this.getNumEdges());
        vertices[i] = i;
      }
      degreeDist = new EnumeratedIntegerDistribution(rng,
          vertices, degreesNormalized);
    }
    return degreeDist;
  }

  private static final Comparator<String> ID_WEIGHT_STRING_COMPARATOR =
      new IdWeightStringComparator();

  /**
   * Load weighted graph into provided 2d arrays.
   *
   * @param lines iterator over lines.
   */
  private static long loadWeightedGraph(SimpleIterator<String> lines, int[][] nbrs, float[][] wts,
                                       float[] wtedOutDegrees, long expectedNumEdges) {
    Optional<String> lineOpt;
    int nodeId = 0;
    long edgesSoFar = 0;
    assert nbrs != null && wts != null;
    while (true) {
      lineOpt = lines.next();
      if (!lineOpt.isPresent()) {
        break;
      } else {
        String line = lineOpt.get();
        if (line.startsWith("%")) {
          continue;
        }
        line = line.trim();
        // Fields are separated by spaces
        if (nodeId >= nbrs.length) {
          throw new RuntimeException("More adjacency lists in input file than specified on "
              + "first line (" + nbrs.length + ")");
        }
        if (line.isEmpty()) { // vertices with no neighbor
          nbrs[nodeId] = new int[0];
          wts[nodeId++] = new float[0];
        } else { // The rest are adjacency lists
          String[] tokens = line.split("\\s+");
          assert tokens.length % 2 == 0;
          int numNeighbors = tokens.length / 2;
          String[] combinedTokensForSorting = new String[numNeighbors];
          for (int i = 0; i < numNeighbors; i++) {
            combinedTokensForSorting[i] = tokens[2 * i] + " " + tokens[2 * i + 1];
          }

          int[] myNbrs = new int[numNeighbors];
          float[] myWts = new float[numNeighbors];
          wtedOutDegrees[nodeId] = 0;
          Arrays.sort(combinedTokensForSorting, ID_WEIGHT_STRING_COMPARATOR);
          for (int i = 0; i < numNeighbors; i++) {
            int nId = getIdFromIdWeightString(combinedTokensForSorting[i]) - 1; // 0-based indexing
            if (nId < 0) {
              throw new IllegalStateException(
                  String.format(
                      "Line %d: neighbor %d is less than 1 (all vertex ids are positive)\n",
                      nodeId + 1, nId + 1
                  )
              );
            }
            myNbrs[i] = nId;
            myWts[i] = getWeightFromIdWeightString(combinedTokensForSorting[i]);
            wtedOutDegrees[nodeId] += myWts[i];
          }
          nbrs[nodeId] = myNbrs;
          wts[nodeId++] = myWts;
          edgesSoFar += numNeighbors;

          if (nodeId > nbrs.length) {
            throw new RuntimeException("Expected number of nodes " + nbrs.length
                + " is inconsistent with number of adjacency lists");
          }
          if (nodeId % 100000 == 0) {
            System.err.print("\rDone loading " + nodeId + " lines.");
          }
        }
      }
    }
    if (edgesSoFar != expectedNumEdges * 2) {
      System.err.println("Expected " + 2 * expectedNumEdges
          + " (counting both directions), but got " + edgesSoFar + " many edges.");
    }
    return edgesSoFar / 2;
  }

  public Graph(SimpleIterator<String> lines) {
    Optional<String> lineOpt;

    lineOpt = lines.next();
    if (!lineOpt.isPresent()) {
      throw new RuntimeException("No lines in input!");
    } else {
      String[] tokens = lineOpt.get().trim().split("\\s+");
      assert tokens.length == 2 || (tokens.length == 3 && tokens[2].equals("1"));
      this.numVertices = Integer.parseInt(tokens[0]);
      long expectedNumEdges = Long.parseLong(tokens[1]);
      this.neighbors = new int[this.numVertices][];
      if (tokens.length > 2 && tokens[2].equals("1")) {
        this.weights = new float[this.numVertices][];
        this.weightedOutDegrees = new float[this.numVertices];
        this.numEdges = loadWeightedGraph(lines,
            this.neighbors, this.weights, this.weightedOutDegrees, expectedNumEdges);
        globalAvgWeight = 0;
        for (float w : this.weightedOutDegrees) {
          globalAvgWeight += w;
        }
        globalAvgWeight = globalAvgWeight / this.numEdges / 2;
      } else {
        this.weights = null;
        this.numEdges = loadUnweightedGraph(lines, this.neighbors, expectedNumEdges);
      }
    }
    checkSizesOfNbrsAndWts();
    checkSymmetryAndSorting();
  }

  /**
   * Returns a Metis representation of the graph, wrapped in an iterable.
   * @param weightFormatter A decimal format that callers can use to specify precision needed
   * in the weights.
   * @return an iterable that returns a freshly constructed iterator each time iterator() is called.
   */
  public Iterable<String> iterableStringRepresentation(DecimalFormat weightFormatter) {
    return new Iterable<String>() {
      // Not thread-safe for now
      @Override
      public Iterator<String> iterator() {
        return new Iterator<String>() {
          private int index = -1;

          public boolean hasNext() {
            return index < numVertices;
          }

          public String next() {
            if (!hasNext()) {
              return null;
            }

            if (index == -1) {
              index++;
              return "" + numVertices + " " + numEdges + (isWeighted() ? " 1" : "");
            }

            StringBuilder sb = new StringBuilder();
            Iterator<WeightedId> wni = getWeightedNeighborsIterator(index);
            boolean isFirst = true;
            while (wni.hasNext()) {
              WeightedId wn = wni.next();
              if (!isFirst) {
                sb.append(" ");
              }
              sb.append(wn.neighborId + 1); // go from 0-indexing to 1-indexing
              if (isWeighted()) {
                sb.append(" ").append(weightFormatter.format(wn.weight));
              }
              isFirst = false;
            }
            index++;
            return sb.toString();

          }
        };
      }
    };
  }

  /**
   * Constructor.
   * N.B. Vertex id's in Metis format is 1-based. We use 0-based indexing.
   *
   * @param filename (required) file in Metis format
   */
  public static Graph fromFile(String filename) throws IOException {
    Graph ret = null;
    BufferedReader br = new BufferedReader(new FileReader(filename));
    SimpleIterator<String> lines = new FileLineIterator(br);
    ret = new Graph(lines);
    br.close();
    return ret;
  }

  /**
   * Get number of vertices.
   *
   * @return number of vertices
   */
  public int getNumVertices() {
    return this.numVertices;
  }

  /**
   * Get all vertices in the graph
   * @return array with ids of all the vertices in the graph
   */
  public int[] getAllVertexIds() {
    if (allVertexIds == null) {
      allVertexIds = new int[this.numVertices];
      for (int i = 0; i < this.numVertices; i++) {
        allVertexIds[i] = i;
      }
    }
    return allVertexIds;
  }

  /**
   * Get number of edges.
   *
   * @return number of edges
   */
  public long getNumEdges() {
    return this.numEdges;
  }

  /**
   * Get neighbors of a vertex.
   *
   * @param i (required) vertex id
   * @return set of neighbors of vertex i
   */
  public int[] getNeighbors(int i) {
    return this.neighbors[i];
  }

  class WeightedNeighborIterator implements Iterator<WeightedId> {
    private int nodeId;
    private int index;
    private int[] myNeighbors;
    private float[] myWeights;

    WeightedNeighborIterator(int nId) {
      nodeId = nId;
      index = 0;
      myNeighbors = getNeighbors(nodeId);
      if (weights != null && weights[nodeId] != null) {
        myWeights = weights[nodeId];
      }
    }

    @Override
    public synchronized boolean hasNext() {
      return index < myNeighbors.length;
    }

    @Override
    public synchronized WeightedId next() {
      float weight;
      int nbr;
      if (myWeights != null) {
        weight = myWeights[index];
      } else {
        weight = 1.0f;
      }
      nbr = myNeighbors[index];
      index++;
      return new WeightedId(nbr, weight);
    }
  }

  /**
   * Get total weight from vertex to a set of vertices
   * @param vertexId source vertex
   * @param set sink vertices
   * @return total weight from vertexId to set
   */
  public double getWeightFromVertexToSet(int vertexId, Set<Integer> set) {
    double ret = 0;
    for (int i : set) {
      ret += getWeightOfEdge(vertexId, i);
    }
    return ret;
  }

  /**
   * Returns array of neighbors for vertex i, sorted in descending order by weight
   */
  public int[] getNeighborsSortedByWeight(int i) {
    if (isWeighted()) {
      int[] indices = new int[neighbors[i].length];
      for (int j = 0; j < indices.length; j++) {
        indices[j] = j;
      }

      int[] sortedIndices =
          Arrays
              .stream(indices)
              .boxed()
              .sorted((id1, id2) -> Float.compare(weights[i][id2], weights[i][id1]))
              .mapToInt(x -> x)
              .toArray();

      int[] ret = new int[neighbors[i].length];
      for (int j = 0; j < indices.length; j++) {
        ret[j] = neighbors[i][sortedIndices[j]];
      }

      return ret;
    } else {
      return getNeighbors(i);
    }
  }

  public Iterator<WeightedId> getWeightedNeighborsIterator(int i) {
    return new WeightedNeighborIterator(i);
  }

  /**
   * Get degree of a vertex.
   *
   * @param i (required) vertex id
   * @return degree of vertex i
   */
  public int getDegree(int i) {
    return this.neighbors[i].length;
  }

  /**
   * Sum of weights on the edges of node i.
   *
   * @param i vertex id
   * @return Sum of weights on the edges of node i.
   */
  public float getWeightedOutDegree(int i) {
    if (this.weightedOutDegrees != null) {
      return this.weightedOutDegrees[i];
    } else {
      return getDegree(i);
    }
  }

  /**
   * Check if two vertices are neighbors by doing binary search over neighbor id's.
   * Search on the smaller array of neighbors.
   *
   * @param i (required) vertex id
   * @param j (required) vertex id
   * @return true iff two vertices share an edge
   */
  public boolean isNeighbors(int i, int j) {
    if (this.getDegree(i) < this.getDegree(j)) {
      return Arrays.binarySearch(this.neighbors[i], j) >= 0;
    } else {
      return Arrays.binarySearch(this.neighbors[j], i) >= 0;
    }
  }

  private float getWeightOfJInI(int i, int j) {
    int indexInI = Arrays.binarySearch(this.neighbors[i], j);
    if (indexInI < 0) {
      return 0;
    } else {
      if (this.weights == null) {
        return 1;
      } else {
        return this.weights[i][indexInI];
      }
    }
  }

  /**
   * Returns weight of the edge (i,i). If no edge exists, then 0 is returned.
   */
  public float getWeightOfEdge(int i, int j) {
    if (this.getDegree(i) < this.getDegree(j)) {
      return getWeightOfJInI(i, j);
    } else {
      return getWeightOfJInI(j, i);
    }
  }

  // Debug
  public void print() {
    System.out.println("numVertices = " + this.numVertices);
    System.out.println("numEdges = " + this.numEdges);
    System.out.println("Graph density: " + this.getDensity());
  }

  /**
   * Creates a Graph object that represents a subgraph of this graph, i.e. which
   * contains only the edges among the specified nodes
   * Note this does not preserve the weights, since there's no usecase for that right now.
   *
   * @param newToOriginal map with mapping from new ids to original ids, the ids need to be from
   * 1 to n
   * @return subgraph among specified nodes
   */
  Graph getUnweightedSubGraph(ImmutableMap<Integer, Integer> newToOriginal) {
    int numNodesOfSubgraph = newToOriginal.size();
    Map<Integer, Integer> originalToNew = new HashMap<>(newToOriginal.size(), 1.0f);

    for (int i = 0; i < numNodesOfSubgraph; i++) {
      originalToNew.put(newToOriginal.get(i), i);
    }

    int[][] newGraphNeighbors = new int[numNodesOfSubgraph][];
    long newGraphEdges = 0;
    for (int i = 0; i < numNodesOfSubgraph; i++) {
      int originalId = newToOriginal.get(i);
      ArrayList<Integer> newNeighbors = new ArrayList<>();
      int[] originalNeighbors = this.getNeighbors(originalId);
      for (int nbr : originalNeighbors) {
        if (originalToNew.containsKey(nbr)) {
          newNeighbors.add(originalToNew.get(nbr));
        }
      }
      newGraphNeighbors[i] = new int[newNeighbors.size()];
      int j = 0;
      for (int newNbr : newNeighbors) {
        newGraphNeighbors[i][j++] = newNbr;
        newGraphEdges++;
      }
    }

    return new Graph(numNodesOfSubgraph, newGraphEdges / 2, newGraphNeighbors);
  }

  /**
   * Score(nodeId, set) = weighted Jaccard similarity between nodeId's neighbors and set.
   * @param nodeId
   * @param set
   * @return
   */
  public double affiliationOfVertexToSet(int nodeId, IntSet set) {
    if (set.size() == 1 && set.contains(nodeId)) {
      return 1.0;
    }

    double wtOfNeighborsInSet = 0;
    int numberOfNeighborsInSet = 0;
    double wtOfNeighborsNotInSet = 0;
    Iterator<WeightedId> iter = this.getWeightedNeighborsIterator(nodeId);
    while (iter.hasNext()) {
      WeightedId wtId = iter.next();
      if (set.contains(wtId.neighborId)) {
        wtOfNeighborsInSet += wtId.weight;
        numberOfNeighborsInSet++;
      } else {
        wtOfNeighborsNotInSet += wtId.weight;
      }
    }
    int setSize = set.contains(nodeId) ? set.size() - 1 : set.size();
    double wtOfSetMinusNeighbors =
        this.getGlobalAvgWeight() * (setSize - numberOfNeighborsInSet);

    double den = wtOfNeighborsInSet + wtOfNeighborsNotInSet + wtOfSetMinusNeighbors;
    if (den == 0.0) {
      return 0.0;
    } else {
      return wtOfNeighborsInSet / den;
    }
  }

  /**
   * Method that retains ONLY the top k edges (in terms of weight) for every vertex in a weighed
   * graph. In case of a conflict (say an edge is in top-k for one of its end points but not for the
   * other) then an OR operation is used, i.e. the edge is retained if it's in the top-k for at
   * least one of its endpoints. Note: method returns new Graph object and does not modify existing
   * object. If the graph is unweighted the method will return the calling object.
   *
   * @param k
   * @return A new Graph object representing the graph after the top-k retention operation.
   */
  public Graph retainTopKEdgesInGraph(int k) {
    //If graph is unweighted, return this object
    if (!this.isWeighted()) {
      return this;
    }
    //Find top-k edges for every vertex using max-heaps, create empty Set for every vertex, and
    //add top-k edges to it
    ArrayList<HashSet<Integer>> topKForVertices = new ArrayList<HashSet<Integer>>();
    for (int vertex = 0; vertex < this.numVertices; vertex++) {
      //Create list of indices of the neighbors in vertex's adjacency list.
      ArrayList<Integer> adjIndexList = new ArrayList<Integer>(this.getDegree(vertex));
      for (int neighborId = 0; neighborId < this.getDegree(vertex); neighborId++) {
        adjIndexList.add(neighborId);
      }
      //sort indices by weights of the edges they correspond to
      float[] wtsOfNeighborsofVertex = this.weights[vertex];
      Collections.sort(adjIndexList,
          (i1, i2) -> (int) Math.signum(wtsOfNeighborsofVertex[i2] - wtsOfNeighborsofVertex[i1]));
      //Create a HashSet containing the top-k edges (neighbors). List contains vertex IDs and not
      //just indices
      HashSet<Integer> topKSet = new HashSet<Integer>();
      int[] vertexNeighbors = this.neighbors[vertex];
      for (int neighbor = 0; neighbor < Math.min(k, this.getDegree(vertex)); neighbor++) {
        topKSet.add(vertexNeighbors[adjIndexList.get(neighbor)]);
      }
      topKForVertices.add(topKSet);
    }
    //Create new adjacency lists and weight lists for creating a new Graph object
    int[][] newAdjLists = new int[this.numVertices][];
    float[][] newWtLists = new float[this.numVertices][];
    //number of edges in new graph. Will initially double count and then divide by two later
    long newNumEdges = 0;
    //Populate using top-k sets computed earlier
    for (int vertex = 0; vertex < this.numVertices; vertex++) {
      //lists for current vertex
      ArrayList<Integer> vertexAdjList = new ArrayList<Integer>();
      ArrayList<Float> vertexWtList = new ArrayList<Float>();
      int[] vertexNeighbors = this.neighbors[vertex];
      for (int oldListIndex = 0; oldListIndex < this.getDegree(vertex); oldListIndex++) {
        int neighborOfVertex = vertexNeighbors[oldListIndex];
        //check if this edge is contained in the top-k of current vertex
        boolean topKinMyList = topKForVertices.get(vertex).contains(neighborOfVertex);
        //check if this edges is contained in the top-k list of the neighbor of the current vertex
        boolean topKinTheirList = topKForVertices.get(neighborOfVertex).contains(vertex);
        if (topKinMyList || topKinTheirList) {
          vertexAdjList.add(neighborOfVertex);
          vertexWtList.add(this.weights[vertex][oldListIndex]);
        }
      }
      //update number of edges
      newNumEdges += vertexAdjList.size();
      //convert ArrayLists to primitive arrays and add to newAdjLists and newWtLists
      newAdjLists[vertex] = new int[vertexAdjList.size()];
      newWtLists[vertex] = new float[vertexAdjList.size()];
      for (int itr = 0; itr < vertexAdjList.size(); itr++) {
        newAdjLists[vertex][itr] = vertexAdjList.get(itr);
        newWtLists[vertex][itr] = vertexWtList.get(itr);
      }
    }
    return new Graph(this.numVertices, newNumEdges / 2, newAdjLists, newWtLists);
  }

  /**
   * Method that squares all the edge weights in the graph. Does not return a new Graph but performs
   * operation in existing Graph object.
   */
  public void squareWeights() {
    //If graph is unweighted, do nothing
    if (this.isWeighted()) {
      //Iterate through weights and square them
      for (int i = 0; i < this.numVertices; i++) {
        for (int j = 0; j < this.getDegree(i); j++) {
          this.weights[i][j] *= this.weights[i][j];
        }
      }
    }
  }
}

