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

import org.junit.Test;

import com.twitter.sbf.core.SparseBinaryMatrix;
import com.twitter.sbf.core.SparseRealMatrix;

import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;

import static org.junit.Assert.assertArrayEquals;

public class SparseRealMatrixTest {
  @Test
  public void normalizeTest() {
    IntSet[] cols = new IntSet[2];
    int[] col0 = {0, 1};
    int[] col1 = {2, 3, 4};
    cols[0] = new IntArraySet(col0);
    cols[1] = new IntArraySet(col1);

    SparseBinaryMatrix z = new SparseBinaryMatrix(5, cols);
    double[][] rowValues = new double[5][];
    for (int i = 0; i < 5; i++) {
      rowValues[i] = new double[1];
    }
    rowValues[0][0] = 1.0;
    rowValues[1][0] = 2.0;
    rowValues[2][0] = Math.sqrt(2);
    rowValues[3][0] = Math.sqrt(0.5);
    rowValues[4][0] = Math.sqrt(0.5);

    SparseRealMatrix srm = new SparseRealMatrix(z, rowValues);
    srm.normalizeToUnitColumn();
    double[] expected0 = {1 / Math.sqrt(5)};
    double[] expected1 = {2 / Math.sqrt(5)};
    assertArrayEquals(expected0, srm.getValuesForRow(0), 1e-5);
    assertArrayEquals(expected1, srm.getValuesForRow(1), 1e-5);

    double[] expected2 = {Math.sqrt(2) / Math.sqrt(3)};
    double[] expected3And4 = {Math.sqrt(0.5) / Math.sqrt(3)};
    assertArrayEquals(expected2, srm.getValuesForRow(2), 1e-5);
    assertArrayEquals(expected3And4, srm.getValuesForRow(3), 1e-5);
    assertArrayEquals(expected3And4, srm.getValuesForRow(4), 1e-5);
  }

  @Test
  public void rowNormalizeTest() {
    IntSet[] cols = new IntSet[2];
    int[] col0 = {0, 1, 2, 4};
    int[] col1 = {2, 3, 4};
    cols[0] = new IntArraySet(col0);
    cols[1] = new IntArraySet(col1);

    SparseBinaryMatrix z = new SparseBinaryMatrix(5, cols);
    double[][] rowValues = new double[5][];
    for (int i = 0; i < 5; i++) {
      if (i == 2 || i == 4) {
        rowValues[i] = new double[2];
      } else {
        rowValues[i] = new double[1];
      }
    }
    rowValues[0][0] = 1.0;
    rowValues[1][0] = 2.0;
    rowValues[2][0] = Math.sqrt(2);
    rowValues[2][1] = Math.sqrt(2);
    rowValues[3][0] = Math.sqrt(0.5);
    rowValues[4][0] = Math.sqrt(1.0);
    rowValues[4][1] = Math.sqrt(1.0);

    SparseRealMatrix srm = new SparseRealMatrix(z, rowValues);
    srm.normalizeToUnitRow();
    double[] expected0 = {1.0};
    double[] expected1 = {1.0};
    assertArrayEquals(expected0, srm.getValuesForRow(0), 1e-5);
    assertArrayEquals(expected1, srm.getValuesForRow(1), 1e-5);

    double[] expected2 = {Math.sqrt(2) / 2.0, Math.sqrt(2) / 2.0};
    double[] expected3 = {1.0};
    double[] expected4 = {1.0 / Math.sqrt(2), 1.0 / Math.sqrt(2)};
    assertArrayEquals(expected2, srm.getValuesForRow(2), 1e-5);
    assertArrayEquals(expected3, srm.getValuesForRow(3), 1e-5);
    assertArrayEquals(expected4, srm.getValuesForRow(4), 1e-5);
  }

  @Test
  public void testApplyCoordinateWiseLogTransform() {
    IntSet[] cols = new IntSet[2];
    int[] col0 = {0, 1};
    int[] col1 = {2, 3, 4};
    cols[0] = new IntArraySet(col0);
    cols[1] = new IntArraySet(col1);

    SparseBinaryMatrix z = new SparseBinaryMatrix(5, cols);
    double[][] rowValues = new double[5][];
    for (int i = 0; i < 5; i++) {
      rowValues[i] = new double[1];
    }
    rowValues[0][0] = 9.0;
    rowValues[1][0] = 99.0;
    rowValues[2][0] = 999.0;
    rowValues[3][0] = 9999.0;
    rowValues[4][0] = 99999.0;

    SparseRealMatrix srm = new SparseRealMatrix(z, rowValues);
    srm.applyCoordinateWiseLogTransform();
    double[] expected0 = {1.0};
    double[] expected1 = {2.0};
    double[] expected2 = {3.0};
    double[] expected3 = {4.0};
    double[] expected4 = {5.0};
    assertArrayEquals(expected0, srm.getValuesForRow(0), 1e-5);
    assertArrayEquals(expected1, srm.getValuesForRow(1), 1e-5);
    assertArrayEquals(expected2, srm.getValuesForRow(2), 1e-5);
    assertArrayEquals(expected3, srm.getValuesForRow(3), 1e-5);
    assertArrayEquals(expected4, srm.getValuesForRow(4), 1e-5);
  }
}
