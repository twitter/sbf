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

/**
 * Class that represents a sparse *non-binary* matrix.
 */
public class SparseRealMatrix {
  // Indicates which (i,j) pairs have non-zero values in this matrix
  private SparseBinaryMatrix nonZeroStructure;

  // realValues[i].length == nonZeroStructure.getRow(i).length
  // realValues[i][j] gives the value for (i, nonZeroStructure.getRow(i)[j]) entry
  private double[][] realValues;
  //The L2-norms of the rows and columns
  //Computed and stored at time of construction of the object.
  private double[] columnNorms;
  private double[] rowNorms;

  public SparseRealMatrix(SparseBinaryMatrix nonZeroStructure, double[][] realValues) {
    this.nonZeroStructure = nonZeroStructure;
    this.realValues = realValues;
    //Initialize column and row norms
    this.columnNorms = new double[nonZeroStructure.getNumCols()];
    this.rowNorms = new double[nonZeroStructure.getNumRows()];
    for (int i = 0; i < nonZeroStructure.getNumRows(); i++) {
      int[] rowWithColIds = nonZeroStructure.getRow(i);
      double[] rowWithValues = realValues[i];
      for (int j = 0; j < rowWithColIds.length; j++) {
        this.columnNorms[rowWithColIds[j]] += rowWithValues[j] * rowWithValues[j];
        this.rowNorms[i] += rowWithValues[j] * rowWithValues[j];
      }
    }
    //Apply square root function to each entry of rowNorms and columnNorms
    for (int row = 0; row < rowNorms.length; row++) {
      rowNorms[row] = Math.sqrt(rowNorms[row]);
    }
    for (int col = 0; col < columnNorms.length; col++) {
      columnNorms[col] = Math.sqrt(columnNorms[col]);
    }
  }

  public int[] getColIdsForRow(int rowId) {
    return nonZeroStructure.getRow(rowId);
  }

  public int getNumRows() {
    return nonZeroStructure.getNumRows();
  }

  public int getNumCols() {
    return nonZeroStructure.getNumCols();
  }

  public double[] getValuesForRow(int rowId) {
    return realValues[rowId];
  }

  public double getRowNorm(int rowId) {
    return rowNorms[rowId];
  }

  public double getColumnNorm(int colId) {
    return columnNorms[colId];
  }

  /**
   * Scale each column in this matrix so that each column has unit length.
   */
  public void normalizeToUnitColumn() {
    for (int i = 0; i < nonZeroStructure.getNumRows(); i++) {
      int[] rowWithColIds = nonZeroStructure.getRow(i);
      for (int j = 0; j < realValues[i].length; j++) {
        realValues[i][j] = realValues[i][j] / columnNorms[rowWithColIds[j]];
      }
    }
  }

  /**
   * Scale each row in this matrix so that each row has unit length.
   */
  public void normalizeToUnitRow() {
    for (int i = 0; i < nonZeroStructure.getNumRows(); i++) {
      int[] rowWithColIds = nonZeroStructure.getRow(i);
      for (int j = 0; j < realValues[i].length; j++) {
        realValues[i][j] = realValues[i][j] / rowNorms[i];
      }
    }
  }

  /**
   * Method to apply a coordinate-wise log10(1+x) transform to the matrix
   */
  public void applyCoordinateWiseLogTransform() {
    for (int i = 0; i < nonZeroStructure.getNumRows(); i++) {
      int[] rowWithColIds = nonZeroStructure.getRow(i);
      for (int j = 0; j < realValues[i].length; j++) {
        realValues[i][j] = Math.log10(1 + realValues[i][j]);
      }
    }
  }

  /**
   * Method to get average nnz per row
   */
  public double getAverageNNZ() {
    return this.nonZeroStructure.nnzPerRow();
  }
}
