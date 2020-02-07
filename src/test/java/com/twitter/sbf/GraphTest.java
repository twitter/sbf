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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.graph.WeightedId;
import com.twitter.sbf.util.SimpleIteratorFromIterator;

import it.unimi.dsi.fastutil.ints.IntArraySet;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class GraphTest {

  public static void testsCommonToWeightedAndUnweighted(Graph g) {
    int[] expectedFor0 = {1, 2, 3, 4};
    assertTrue(Arrays.equals(g.getNeighbors(0), expectedFor0));

    int[] expectedFor7 = {4, 5, 6};
    assertTrue(Arrays.equals(g.getNeighbors(7), expectedFor7));

    assertTrue(g.getNumEdges() == 13L);
    assertTrue(g.getNumVertices() == 9);
    assertTrue(g.isNeighbors(1, 2));
    assertFalse(g.isNeighbors(1, 7));
    assertFalse(g.isNeighbors(2, 6));
    assertTrue(g.getDegree(4) == 4);
    assertTrue(g.getNeighbors(8).length == 0);

    assertTrue(Math.abs(g.getDensity() - 2 * 13.0f / 9 / 8) < 1e-5);
  }

  private static Graph loadUnweighted() {
    List<String> lines = new ArrayList<>();
    lines.add("9 13");
    lines.add("2 3 4 5"); // 0
    lines.add("1 3 4"); // 1
    lines.add("1 2 4"); // 2
    lines.add("1 2 3"); // 3
    lines.add("1 6 7 8"); // 4
    lines.add("7 5 8"); // 5
    lines.add("5 6 8"); // 6
    lines.add("7 6 5"); // 7
    lines.add(""); // 8 is a singleton node

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  private static Graph loadSimpleWeighted() {
    List<String> lines = new ArrayList<>();
    lines.add("8 12 1");
    lines.add("2 0.5 3 0.5 4 0.5 5 1.0"); // 0
    lines.add("1 0.5 3 0.5 4 0.5"); // 1
    lines.add("1 0.5 2 0.5 4 0.5"); // 2
    lines.add("1 0.5 2 0.5 3 0.5"); // 3
    lines.add("1 1.0 6 1.0 7 1.0 8 1.0"); // 4
    lines.add("5 1.0 7 1.0"); // 5
    lines.add("5 1.0 6 1.0 8 1.0"); // 6
    lines.add("5 1.0 7 1.0"); // 7

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  private static Graph loadWeighted() {
    List<String> lines = new ArrayList<>();
    lines.add("9 13 1");
    lines.add("2 0.2 3 0.3 4 5e-2 5 1e-3"); // 0
    lines.add("1 0.2 3 0.3 4 0.4"); // 1
    lines.add("1 0.3 2 0.3 4 0.25"); // 2
    lines.add("1 0.05 2 0.4 3 0.25"); // 3
    lines.add("1 0.001 6 0.15 7 0.3 8 0.2"); // 4
    lines.add("5 0.15 7 0.2 8 0.3"); // 5
    lines.add("5 0.3 6 0.2 8 0.01"); // 6
    lines.add("5 0.2 6 0.3 7 1e-2"); // 7
    lines.add(""); // 8 is a singleton node

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  @Test
  public void testLoadingUnweighted() {
    Graph g = loadUnweighted();
    testsCommonToWeightedAndUnweighted(g);

    assertTrue(Math.abs(g.getWeightedOutDegree(4) - 4.0f) < 1e-5);
    assertTrue(Math.abs(g.getWeightedOutDegree(8) - 0f) < 1e-5);

    Iterator<WeightedId> iterFor4 = g.getWeightedNeighborsIterator(5);
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(4, 1f)));
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(6, 1f)));
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(7, 1f)));
  }

  @Test
  public void testLoadingWeighted() {
    Graph g = loadWeighted();
    testsCommonToWeightedAndUnweighted(g);
    assertTrue(Math.abs(g.getWeightedOutDegree(0) - 0.551) < 1e-5);
    assertTrue(Math.abs(g.getWeightedOutDegree(7) - 0.51) < 1e-5);

    Iterator<WeightedId> iterFor4 = g.getWeightedNeighborsIterator(4);
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(0, 0.001f)));
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(5, 0.15f)));
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(6, 0.3f)));
    assertTrue(iterFor4.next().equalsWithTolerance(new WeightedId(7, 0.2f)));
  }

  @Test
  public void testConductance() {
    Graph g = loadUnweighted();

    double c1 = g.getNeighborhoodConductanceUnweighted(1);
    double c4 = g.getNeighborhoodConductanceUnweighted(4);
    assertTrue(Math.abs(c1 - 1.0 / 13.0) < 1e-5);
    assertTrue(Math.abs(c4 - 3.0 / 9.0) < 1e-5);
    assertTrue(Math.abs(g.getNeighborhoodConductance(4) - 3.0 / 9.0) < 1e-5);
    assertTrue(Math.abs(g.getNeighborhoodConductance(1) - 1.0 / 13.0) < 1e-5);

    Graph gWted = loadWeighted();

    double c1Weighted = gWted.getNeighborhoodConductance(1);
    assertTrue(Math.abs(c1Weighted - 1e-3 / 2.321) < 1e-5);
  }

  @Test
  public void testGetNeighborsSortedByWeight() {
    Graph g = loadWeighted();
    int[] expected = {3, 2, 0};
    assertArrayEquals(g.getNeighborsSortedByWeight(1), expected);

    int[] expected2 = {6, 7, 5, 0};
    assertArrayEquals(g.getNeighborsSortedByWeight(4), expected2);
  }

  @Test
  public void testAffiliationOfVerticesToSets() {
    Graph g = loadSimpleWeighted();
    assertTrue(Math.abs(g.getGlobalAvgWeight() - 0.75) < 1e-5);

    int[] set = {0, 1, 2, 3, 7};
    double aff = g.affiliationOfVertexToSet(0, new IntArraySet(set));
    assertTrue(Math.abs(aff - 1.5 / 3.25) < 1e-5);
  }

  @Test
  public void testRetainTopKEdgesInGraph() {
    //Test with simple weighted graph
    Graph g = loadSimpleWeighted();
    Graph gTopK = g.retainTopKEdgesInGraph(2);
    assertTrue(Math.abs(gTopK.getWeightOfEdge(2,  3) - 0.0) < 1e-5);
    assertTrue(gTopK.getDegree(7) == 2);
    assertTrue(gTopK.getDegree(3) == 2);
    assertTrue(gTopK.getDegree(4) == 4);
    assertTrue(gTopK.getNumEdges() == 11);
    //Test with weighted graph;
    Graph gHard = loadWeighted();
    Graph gHardTopK = gHard.retainTopKEdgesInGraph(2);
    //expected neighbors of vertex 1, 3, 5, 8; 0-indexed inside in the array
    int[] expected1 = {1, 2};
    int[] expected3 = {0, 1, 3};
    int[] expected5 = {6, 7};
    int[] expected8 = {4, 5};
    assertArrayEquals(expected1, gHardTopK.getNeighbors(0));
    assertArrayEquals(expected3, gHardTopK.getNeighbors(2));
    assertArrayEquals(expected5, gHardTopK.getNeighbors(4));
    assertArrayEquals(expected8, gHardTopK.getNeighbors(7));
  }

  @Test
  public void testSquareWeightsInGraph() {
    Graph g = loadWeighted();
    g.squareWeights();
    assertTrue(Math.abs(g.getWeightOfEdge(0, 4) - 1e-6) < 1e-10);
    assertTrue(Math.abs(g.getWeightOfEdge(6, 4) - 0.09) < 1e-5);
    assertTrue(Math.abs(g.getWeightOfEdge(2, 3) - 0.0625) < 1e-5);
  }
}

