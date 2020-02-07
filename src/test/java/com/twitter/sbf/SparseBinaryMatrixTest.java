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

import java.util.List;
import java.util.Set;

import com.google.common.collect.ImmutableList;

import org.junit.Test;

import com.twitter.sbf.core.SparseBinaryMatrix;
import com.twitter.sbf.util.SimpleIteratorFromIterator;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

public class SparseBinaryMatrixTest {

  private boolean compareSetAndArray(Set<Integer> a, int[] arr) {
    return a.containsAll(new IntOpenHashSet(arr)) && a.size() == arr.length;
  }

  @Test
  public void testInitFromRows() {
    SparseBinaryMatrix z = new SparseBinaryMatrix(5, 2);

    List<String> rows = ImmutableList.of(
      "1 4",
      "",
      "2",
      "1 2 4",
      "1"
    );

    z.initFromRows(new SimpleIteratorFromIterator<>(rows.iterator()));

    int[] col1 = {0, 3, 4};
    int[] col2 = {2, 3};
    int[] col3 = new int[0];
    int[] col4 = {0, 3};

    assertTrue(compareSetAndArray(z.getColumn(0), col1));
    assertTrue(compareSetAndArray(z.getColumn(1), col2));
    assertTrue(compareSetAndArray(z.getColumn(2), col3));
    assertTrue(compareSetAndArray(z.getColumn(3), col4));
    assertTrue(z.getNumCols() == 4);

    int[] row1 = {0, 3};
    int[] row2 = {};
    int[] row3 = {1};
    int[] row4 = {0, 1, 3};
    int[] row5 = {0};
    assertArrayEquals(z.getRow(0), row1);
    assertArrayEquals(z.getRow(1), row2);
    assertArrayEquals(z.getRow(2), row3);
    assertArrayEquals(z.getRow(3), row4);
    assertArrayEquals(z.getRow(4), row5);

  }
}
