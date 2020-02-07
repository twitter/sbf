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
import java.util.Arrays;
import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomAdaptor;
import org.junit.Test;

import com.twitter.sbf.core.AlgorithmConfig;
import com.twitter.sbf.core.IdCountWtContributors;
import com.twitter.sbf.core.IdCountWtTriple;
import com.twitter.sbf.core.MHAlgorithm;
import com.twitter.sbf.core.SparseBinaryMatrix;
import com.twitter.sbf.core.SparseRealMatrix;
import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.util.SimpleIteratorFromIterator;

import it.unimi.dsi.fastutil.ints.IntSet;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

public class AlgorithmTest {
  private static Graph simple8NodeGraph() {
    List<String> lines = new ArrayList<>();
    lines.add("8 13");
    lines.add("2 3 4 5"); // 0
    lines.add("1 3 4"); // 1
    lines.add("1 2 4"); // 2
    lines.add("1 2 3"); // 3
    lines.add("1 6 7 8"); // 4
    lines.add("5 7 8"); // 5
    lines.add("5 6 8"); // 6
    lines.add("5 6 7"); // 7

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  private static Graph simpleWeighted8NodeGraph() {
    List<String> lines = new ArrayList<>();
    lines.add("8 13 1");
    lines.add("2 1.0 3 1.0 4 1.0 5 7.0"); // 0
    lines.add("1 1.0 3 1.0 4 1.0"); // 1
    lines.add("1 1.0 2 1.0 4 1.0"); // 2
    lines.add("1 1.0 2 1.0 3 1.0"); // 3
    lines.add("1 7.0 6 1.0 7 1.0 8 1.0"); // 4
    lines.add("5 1.0 7 1.0 8 1.0"); // 5
    lines.add("5 1.0 6 1.0 8 1.0"); // 6
    lines.add("5 1.0 6 1.0 7 1.0"); // 7

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  private static SparseBinaryMatrix initZForSimple8NodeGraph(Graph g, int k) {
    assertTrue(g.getNumVertices() == 8);
    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), k);
    List<String> initZ = new ArrayList<>();
    initZ.add("2");
    initZ.add("1");
    initZ.add("1");
    initZ.add("1");
    initZ.add("1");
    initZ.add("2");
    initZ.add("2");
    initZ.add("2");

    z.initFromRows(new SimpleIteratorFromIterator<String>(initZ.iterator()));
    return z;
  }

  private void checkResultForSimple8NodeGraph(SparseBinaryMatrix s) {
    for (int i = 0; i < 8; i++) {
      System.out.println("row " + i + ": " + Arrays.toString(s.getRow(i)));
    }

    assertTrue(s.getRow(1).length == 1);
    assertTrue(s.getRow(5).length == 1);

    assertTrue(
        (s.getRow(0).length == 1 && s.getRow(0)[0] == s.getRow(1)[0])
        || (s.getRow(0).length == 2 && (s.getRow(0)[0] == s.getRow(1)[0]
            || s.getRow(0)[1] == s.getRow(1)[0]))
    );
    assertArrayEquals(s.getRow(1), s.getRow(2));
    assertArrayEquals(s.getRow(1), s.getRow(3));

    assertTrue(
        (s.getRow(4).length == 1 && s.getRow(4)[0] == s.getRow(5)[0])
            || (s.getRow(4).length == 2 && (s.getRow(4)[0] == s.getRow(5)[0]
                || s.getRow(4)[1] == s.getRow(5)[0]))
    );
    assertArrayEquals(s.getRow(5), s.getRow(6));
    assertArrayEquals(s.getRow(5), s.getRow(7));
  }

  @Test
  public void testSimpleWithRandomInit() {
    Graph g = simple8NodeGraph();

    AlgorithmConfig cfg =
        new AlgorithmConfig().withK(2).withScaleCoeff(1.0)
            .withRunAllEpochs(true).withMaxEpoch(10).withCpu(2).withEvalRatio(1.0)
            .withRng(new RandomAdaptor(new JDKRandomGenerator(1)));
    cfg.updateImmediately = true;
    SparseBinaryMatrix z = new SparseBinaryMatrix(g.getNumVertices(), cfg.k);
    z.initEmptyRowsRandomly(cfg.rng);
    MHAlgorithm algo = new MHAlgorithm(cfg, g, z, new PrintWriter(System.err));
    checkResultForSimple8NodeGraph(algo.optimize());
  }

  @Test
  public void testSimpleWithNeighborhoodInit() {
    Graph g = simple8NodeGraph();

    AlgorithmConfig cfg =
        new AlgorithmConfig().withK(2).withScaleCoeff(1.0)
            .withRunAllEpochs(true).withMaxEpoch(5).withCpu(1).withEvalRatio(1.0);

    SparseBinaryMatrix sIn = new SparseBinaryMatrix(g.getNumVertices(), cfg.k);
    PrintWriter err = new PrintWriter(System.err);
    sIn.initFromBestNeighborhoods(g, Graph::getNeighborhoodConductance, false, err);
    MHAlgorithm algo = new MHAlgorithm(cfg, g, sIn, err);
    checkResultForSimple8NodeGraph(algo.optimize());

  }

  @Test
  public void testSimpleWithHandPickedInit() {
    Graph g = simple8NodeGraph();
    SparseBinaryMatrix sIn = initZForSimple8NodeGraph(g, 2);
    AlgorithmConfig cfg =
        new AlgorithmConfig().withK(2).withScaleCoeff(1.0)
            .withRunAllEpochs(true).withMaxEpoch(5).withCpu(1).withEvalRatio(1.0);
    MHAlgorithm algo = new MHAlgorithm(cfg, g, sIn, new PrintWriter(System.err));
    checkResultForSimple8NodeGraph(algo.optimize());
  }

  @Test
  public void testSumRows() {
    SparseBinaryMatrix z = new SparseBinaryMatrix(5, 5);
    z.initFromRows(
        new SimpleIteratorFromIterator<>(
            ImmutableList.of("1", "2 5", "3", "4", "5").iterator()
        )
    );
    // just a graph in the shape of a pentagon
    Graph g = new Graph(
        new SimpleIteratorFromIterator<>(
            ImmutableList.of("5 5 1", "2 0.1 5 1.0", "1 0.1 3 0.5",
                "2 0.5 4 0.5", "3 0.5 5 1.0", "1 1.0 4 1.0")
                .iterator()
        )
    );

    IdCountWtTriple[] sumFor0 = MHAlgorithm.sumRows(g.getNeighbors(0), z,
        neighborId -> g.getWeightOfEdge(0, neighborId));
    assertTrue(sumFor0.length == 2);
    assertTrue(sumFor0[0].id == 1);
    assertTrue(sumFor0[0].count == 1);
    assertTrue(Math.abs(sumFor0[0].wt - 0.1) < 1e-5);
    assertTrue(sumFor0[1].id == 4);
    assertTrue(sumFor0[1].count == 2);
    assertTrue(Math.abs(sumFor0[1].wt - 1.1) < 1e-5);

    IdCountWtContributors[] fullSumFor0 = MHAlgorithm.sumRowsWithContributors(
        g.getNeighbors(0), z,
        neighborId -> g.getWeightOfEdge(0, neighborId));
    assertTrue(fullSumFor0.length == 2);
    assertTrue(fullSumFor0[0].id == 1);
    assertTrue(fullSumFor0[0].count == 1);
    assertTrue(fullSumFor0[0].contributors.equals(ImmutableSet.of(1)));
    assertTrue(Math.abs(fullSumFor0[0].wt - 0.1) < 1e-5);
    assertTrue(fullSumFor0[1].id == 4);
    assertTrue(fullSumFor0[1].count == 2);
    assertTrue(Math.abs(fullSumFor0[1].wt - 1.1) < 1e-5);
    assertTrue(fullSumFor0[1].contributors.equals(ImmutableSet.of(1, 4)));
  }

  @Test
  public void testBestPrefixConductanceSet() {
    Graph g = simple8NodeGraph();
    int[] testSet = {0, 1, 2, 3, 4, 5, 6};
    int[] s1 = MHAlgorithm.getBestPrefixConductanceSet(g, testSet);
    int[] expected = {0, 1, 2, 3};
    assertArrayEquals(expected, s1);

    int[] testSet2 = {0, 1, 2};
    int[] s2 = MHAlgorithm.getBestPrefixConductanceSet(g, testSet2);
    int[] expected2 = {0, 1, 2};
    assertArrayEquals(expected2, s2);

    int[] testSet3 = {0, 1, 2, 3, 4, 5};
    int[] s3 = MHAlgorithm.getBestPrefixConductanceSet(simpleWeighted8NodeGraph(), testSet3);
    int[] expected3 = {0, 1, 2, 3, 4};
    assertArrayEquals(expected3, s3);

    int[] testSet4 = {0, 1};
    IntSet[] s = MHAlgorithm.getBestConductanceSubNeighborhoods(g, testSet4);
    int[] expected4 = {1, 2, 3};
    int[] s4 = Arrays.stream(s[0].toArray()).mapToInt(i -> (Integer) i).toArray();
    Arrays.sort(s4);
    Arrays.sort(expected4);
    assertArrayEquals(s4, expected4);
    int[] expected5 = {0, 2, 3};
    int[] s5 = Arrays.stream(s[1].toArray()).mapToInt(i -> (Integer) i).toArray();
    Arrays.sort(s5);
    Arrays.sort(expected5);
    assertArrayEquals(s5, expected5);
  }

  @Test
  public void testClusterAssignmentScores() {
    Graph g = simple8NodeGraph();
    SparseBinaryMatrix z = new SparseBinaryMatrix(8, 2);
    z.initFromRows(
        new SimpleIteratorFromIterator<>(
            ImmutableList.of("1", "1", "1", "1", "2", "2", "2", "2").iterator()
        )
    );
    SparseRealMatrix srm = MHAlgorithm.heuristicallyScoreClusterAssignments(g, z);
    assertTrue(srm.getNumRows() == 8);
    assertArrayEquals(srm.getColIdsForRow(0), srm.getColIdsForRow(1));
    assertArrayEquals(srm.getColIdsForRow(0), srm.getColIdsForRow(2));
    assertArrayEquals(srm.getColIdsForRow(0), srm.getColIdsForRow(3));
    assertArrayEquals(srm.getColIdsForRow(4), srm.getColIdsForRow(5));
    assertArrayEquals(srm.getColIdsForRow(4), srm.getColIdsForRow(6));
    assertArrayEquals(srm.getColIdsForRow(4), srm.getColIdsForRow(7));

    double[] expectedFor6OutOf8 = {1.0};
    double[] expectedFor0And4 = {0.75};
    assertArrayEquals(srm.getValuesForRow(0), expectedFor0And4, 1e-5);
    assertArrayEquals(srm.getValuesForRow(4), expectedFor0And4, 1e-5);
    assertArrayEquals(srm.getValuesForRow(1), expectedFor6OutOf8, 1e-5);
    assertArrayEquals(srm.getValuesForRow(2), expectedFor6OutOf8, 1e-5);
    assertArrayEquals(srm.getValuesForRow(3), expectedFor6OutOf8, 1e-5);
    assertArrayEquals(srm.getValuesForRow(5), expectedFor6OutOf8, 1e-5);
    assertArrayEquals(srm.getValuesForRow(6), expectedFor6OutOf8, 1e-5);
    assertArrayEquals(srm.getValuesForRow(7), expectedFor6OutOf8, 1e-5);
  }

}
