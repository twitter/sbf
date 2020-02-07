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
import java.util.List;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.twitter.sbf.graph.BipartiteGraph;
import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.util.SimpleIteratorFromIterator;

public class BipartiteGraphTest {

  private static BipartiteGraph loadBipartiteGraph() {
    //Bipartite graph with 7 vertices on left, 4 vertices on the right, and 16 edges in total
    //Note the METIS format
    List<String> lines = new ArrayList<>();
    lines.add("7 4 18");
    //Adjacency lists for left side vertices
    lines.add("1 2 4"); // 0
    lines.add("1 2 4"); // 1
    lines.add("1 2"); // 2
    lines.add("2 3 4"); // 3
    lines.add("1 3 4"); // 4
    lines.add("3 4"); // 5
    lines.add("3 4"); // 6
    //Adjacency lists for right side vertices
    lines.add("1 2 3 5"); // 0
    lines.add("1 2 3 4"); // 1
    lines.add("4 5 6 7"); // 2
    lines.add("1 2 4 5 6 7"); //3

    return new BipartiteGraph(new SimpleIteratorFromIterator<>(lines.iterator()));
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
   * Tests the following methods:
   * getNumLeftVertices, getNumRightVertices,
   * getNeighborsForLeftById, getNeighborsForRightByRightId
   * getNumEdges, and getEdgeDensity
   * Implictly tests the constructor of the BipartiteGraph class that reads from SimpleIterator
   */
  @Test
  public void testAccessMethods() {
    BipartiteGraph bg = loadBipartiteGraph();
    //Test edge count
    assertTrue(bg.getNumEdges() == 18);

    //Test number of left vertices, right vertices, and edge density
    assertTrue(bg.getNumLeftVertices() == 7);
    assertTrue(bg.getNumRightVertices() == 4);
    assertTrue(Math.abs(bg.getEdgeDensity() - 18.0f / 7 / 4) < 1e-5);

    //Test getNeighborsForLeftById
    int[] expectedFor2L = {0, 1, 3};
    assertArrayEquals(expectedFor2L, bg.getNeighborsForLeftById(1));
    int[] expectedFor6L = {2, 3};
    assertArrayEquals(expectedFor6L, bg.getNeighborsForLeftById(5));
    int[] expectedFor4L = {1, 2, 3};
    assertArrayEquals(expectedFor4L, bg.getNeighborsForLeftById(3));
    int[] expectedFor5L = {0, 2, 3};
    assertArrayEquals(expectedFor5L, bg.getNeighborsForLeftById(4));

    //Test getNeighborsForRightById
    int[] expectedFor2R = {0, 1, 2, 3};
    assertArrayEquals(expectedFor2R, bg.getNeighborsForRightById(1));
    int[] expectedFor1R = {0, 1, 2, 4};
    assertArrayEquals(expectedFor1R, bg.getNeighborsForRightById(0));
    int[] expectedFor3R = {3, 4, 5, 6};
    assertArrayEquals(expectedFor3R, bg.getNeighborsForRightById(2));
    int[] expectedFor4R = {0, 1, 3, 4, 5, 6};
    assertArrayEquals(expectedFor4R, bg.getNeighborsForRightById(3));
  }

  //Tests the getTranspose method of the BipartiteGraph class
  @Test
  public void testTranspose() {
    BipartiteGraph bg = loadBipartiteGraph();
    BipartiteGraph bgTrans = bg.getTranspose();

    //Check if left side of bg is the right side of bgTrans and vice-versa
    assertFalse(bg.getNumLeftVertices() == bgTrans.getNumLeftVertices());
    assertTrue(bg.getNumLeftVertices() == bgTrans.getNumRightVertices());
    assertTrue(bg.getNumRightVertices() == bgTrans.getNumLeftVertices());

    assertArrayEquals(bg.getNeighborsForLeftById(0), bgTrans.getNeighborsForRightById(0));
    assertArrayEquals(bg.getNeighborsForLeftById(2), bgTrans.getNeighborsForRightById(2));
    assertArrayEquals(bg.getNeighborsForLeftById(6), bgTrans.getNeighborsForRightById(6));

    assertArrayEquals(bg.getNeighborsForRightById(0), bgTrans.getNeighborsForLeftById(0));
    assertArrayEquals(bg.getNeighborsForRightById(3), bgTrans.getNeighborsForLeftById(3));
    assertArrayEquals(bg.getNeighborsForRightById(2), bgTrans.getNeighborsForLeftById(2));
  }

  //Tests if projectRight works
  @Test
  public void testProjectRight() {
    BipartiteGraph bg = loadBipartiteGraph();
    //Test right projection

    Graph bgProjRight = bg.projectRight(0.0);
    //Test for right number of edges and vertices
    assertTrue(bgProjRight.getNumVertices() == 4);
    assertTrue(bgProjRight.getNumEdges() == 6);
    assertTrue(bgProjRight.isWeighted());

    //Check if the edge weights in the projected graph match to expected values
    float expectedFor12 = (float) (3.0f / Math.sqrt(16));
    float expectedFor13 = (float) (1.0f / Math.sqrt(16));
    float expectedFor14 = (float) (3.0f / Math.sqrt(24));
    float expectedFor23 = (float) (1.0f / Math.sqrt(16));
    float expectedFor24 = (float) (3.0f / Math.sqrt(24));
    float expectedFor34 = (float) (4.0f / Math.sqrt(24));
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 0) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 1) - expectedFor12) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 2) - expectedFor13) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 3) - expectedFor14) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(1, 1) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(1, 2) - expectedFor23) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(1, 3) - expectedFor24) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(2, 2) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(2, 3) - expectedFor34) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(3, 3) - 0.0) < 1e-5;
  }


  //Tests it projectLeft works as expected
  @Test
  public void testProjectLeft() {
    BipartiteGraph bg = loadBipartiteGraph();
    Graph bgProjLeft = bg.projectLeft(0.0);

    //Test for right number of edges and vertices
    assertTrue(bgProjLeft.getNumVertices() == 7);
    assertTrue(bgProjLeft.getNumEdges() == 19);

    //Test for edge weights and adjacency
    //vertex 3 on left side does not share any neighbors with vertices 6 and 7 on the left
    //Thus, 3 is not adjacent to 6 and 7 in the leftProjected graph.
    assertFalse(bgProjLeft.isNeighbors(2, 5));
    assertFalse(bgProjLeft.isNeighbors(2, 6));
    //No self loop test
    assertFalse(bgProjLeft.isNeighbors(3, 3));

    //checking weights for some edges
    assertTrue(Math.abs(bgProjLeft.getWeightOfEdge(0, 2) - 2.0f / Math.sqrt(6)) < 1e-5);
    assertTrue(Math.abs(bgProjLeft.getWeightOfEdge(1, 4) - 2.0f / Math.sqrt(9)) < 1e-5);
    assertTrue(Math.abs(bgProjLeft.getWeightOfEdge(1, 6) - 1.0f / Math.sqrt(6)) < 1e-5);
    assertTrue(Math.abs(bgProjLeft.getWeightOfEdge(5, 6) - 1.0f) < 1e-5);
    assertTrue(Math.abs(bgProjLeft.getWeightOfEdge(2, 3) - 1.0f / Math.sqrt(6)) < 1e-5);
  }

  //Tests if projectRight with thresholding works, threshold = 0.7
  @Test
  public void testProjectRightWithThreshold() {
    BipartiteGraph bg = loadBipartiteGraph();
    //Test right projection

    Graph bgProjRight = bg.projectRight(0.7);
    //Test for right number of edges and vertices
    assertTrue(bgProjRight.getNumVertices() == 4);
    assertTrue(bgProjRight.getNumEdges() == 2);

    //Check if the edge weights in the projected graph with threshold = 0.7 match to expected values
    float expectedFor12 = (float) (3.0f / Math.sqrt(16));
    float expectedFor34 = (float) (4.0f / Math.sqrt(24));
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 0) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 1) - expectedFor12) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 2) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(0, 3) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(1, 1) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(1, 2) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(1, 3) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(2, 2) - 0.0) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(2, 3) - expectedFor34) < 1e-5;
    assert Math.abs(bgProjRight.getWeightOfEdge(3, 3) - 0.0) < 1e-5;
  }
}
