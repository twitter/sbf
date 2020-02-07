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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Optional;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.BiFunction;

import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.util.FileLineIterator;
import com.twitter.sbf.util.SimpleIterator;
import com.twitter.sbf.util.Util;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;


/**
 * Sparse numRows x numCols binary matrix, stores latent factors.
 * Rows represent vertices, columns represent communities.
 * Each row is a sorted list for efficient dot product.
 * Each column is a set for efficient union operation.
 * <p>
 */
public class SparseBinaryMatrix {
  private int numRows;
  private int numCols;
  private int[][] rows; // each row is sorted in ascending order
  private IntSet[] cols;

  /**
   * Constructor. Allocates empty rows and columns.
   *
   * @param numRows (required) number of rows
   * @param numCols (required) number of columns
   */
  public SparseBinaryMatrix(int numRows, int numCols) {
    this.numRows = numRows;
    this.numCols = numCols;
    this.rows = new int[this.numRows][];
    for (int i = 0; i < this.numRows; i++) {
      this.rows[i] = new int[0]; //Util.EMPTY_INT_ARRAY;
    }
    this.cols = new IntSet[this.numCols];
  }

  public SparseBinaryMatrix(int numRows, IntSet[] cols) {
    this(numRows, cols.length);
    initFromColSets(cols);
  }

  public int getNumCols() {
    return this.numCols;
  }

  public int getNumRows() {
    return this.numRows;
  }

  /**
   * Reset all rows and columns as empty.
   */
  public void reset() {
    for (int i = 0; i < this.numRows; i++) {
      this.rows[i] = new int[0];
    }
    for (int k = 0; k < this.numCols; k++) {
      this.cols[k] = new IntOpenHashSet();
    }
  }

  /**
   * Set one nonzero randomly for each non-empty row.
   */
  public void initEmptyRowsRandomly(Random rng) {
    this.reset();
    for (int i = 0; i < this.numRows; i++) {
      if (this.rows[i] == null || this.rows[i].length == 0) {
        int randomColId = rng.nextInt(this.numCols);
        this.rows[i] = new int[]{randomColId};
        this.getColumn(randomColId).add(i);
      }
    }
  }

  /**
   * initialize matrix from rows.
   * each row is sparse binary, and has white-space separated indices of columns,
   * between 1 and this.numCols
   */
  public void initFromRows(SimpleIterator<String> lines) {
    this.reset();
    int rowId = 0;
    Optional<String> lineOpt;
    Int2ObjectMap<IntArrayList> colsMap = new Int2ObjectOpenHashMap<>();
    int maxColId = -1;
    while (true) {
      lineOpt = lines.next();
      if (!lineOpt.isPresent()) {
        break;
      } else {
        String line = lineOpt.get().trim();
        if (!line.isEmpty()) {
          String[] tokens = line.split("\\s+");
          int[] row = new int[tokens.length];
          for (int i = 0; i < tokens.length; i++) {
            int colId = Integer.parseInt(tokens[i]) - 1; // convert to 0-based
            if (colId < 0) {
              throw new RuntimeException(
                String.format(
                  "Column %d smaller than 1, in line number %d, line:'%s'",
                    colId + 1, rowId + 1, line
                )
              );
            }
            row[i] = colId;
            if (!colsMap.containsKey(colId)) {
              colsMap.put(colId, new IntArrayList());
            }
            colsMap.get(colId).add(rowId);
            if (colId > maxColId) {
              maxColId = colId;
            }
          }
          Arrays.sort(row);
          this.rows[rowId] = row;
        }
        rowId++;
        if (rowId > this.numRows) {
          throw new RuntimeException(
              "More rows in input rows file than the expected " + this.numRows
          );
        }
      }
    }
    if (maxColId > this.numCols) {
      this.numCols = maxColId + 1;
      this.cols = new IntSet[this.numCols];
    }
    for (int colId = 0; colId < this.numCols; colId++) {
      if (colsMap.containsKey(colId) && !colsMap.get(colId).isEmpty()) {
        this.cols[colId] = new IntOpenHashSet(colsMap.get(colId));
      } else {
        this.cols[colId] = new IntOpenHashSet();
      }
    }
  }

  /**
   * Init from file that contains row view of the matrix.
   * Each line contains nonzero column id's of a row, separated by spaces.
   * Column id's are 1-based.
   *
   * @param filename (required) file where each line is a row
   */
  public void initFromRows(String filename) throws IOException {
    long numLines = Util.getNumLinesOfFile(filename);
    if (numLines != this.numRows) {
      throw new IllegalStateException(
          String.format("Rows file has %d lines, expecting %d", numLines, this.numRows)
      );
    }

    // Read from file
    BufferedReader br = new BufferedReader(new FileReader(filename));
    initFromRows(new FileLineIterator(br));
    br.close();
  }

  /**
   * Initialize matrix from columns, each column represented as
   * the set of row ids for which this column is active.
   */
  public void initFromColSets(IntSet[] inputCols) {
    // Initialize resizable rows
    IntList[] dynamicRows = new IntList[this.numRows];
    for (int i = 0; i < this.numRows; i++) {
      dynamicRows[i] = new IntArrayList();
    }

    this.cols = inputCols;
    for (int i = 0; i < this.numCols; i++) {
      for (int c : this.cols[i]) {
        dynamicRows[c].add(i);
      }
    }

    // Sort each row and copy to myself
    for (int i = 0; i < this.numRows; i++) {
      int[] row = dynamicRows[i].toIntArray();
      Arrays.sort(row);
      this.rows[i] = row;
    }
  }

  /**
   * Initialize the columns for this matrix using the neighborhoods with
   * the least cost function values.
   */
  public void initFromBestNeighborhoods(
      Graph g,
      BiFunction<Graph, Integer, Double> costForNeighborhoodFn,
      boolean allowOverlap,
      PrintWriter diagnosticsWriter
  ) {
    initFromColSets(
        MHAlgorithm.getBottomKNeighborhoods(
            g, this.numCols, costForNeighborhoodFn, allowOverlap, diagnosticsWriter
        )
    );
  }

  /**
   * Write each row as a line in a file.
   * Items are separated by a space, sorted in ascending order.
   * Add 1 back to make it 1-indexed.
   *
   * @param filename (required) output file
   */
  public void outputByRows(String filename) throws IOException {
    BufferedWriter rowsWriter = new BufferedWriter(new FileWriter(new File(filename)));
    for (int i = 0; i < this.numRows; i++) {
      for (int j = 0; j < this.rows[i].length; j++) {
        if (j > 0) {
          rowsWriter.write(" ");
        }
        // Add 1 to get back 1-indexing
        rowsWriter.write(Integer.toString(this.rows[i][j] + 1));
      }
      rowsWriter.newLine();
    }
    rowsWriter.close();
  }

  /**
   * Write each column as a line in a file.
   * Items are separated by a space, unsorted.
   * Add 1 back to make it 1-indexed.
   *
   * @param filename (required) output file
   */
  public void outputByCols(String filename) throws IOException {
    BufferedWriter colsWriter = new BufferedWriter(new FileWriter(new File(filename)));
    for (int k = 0; k < this.numCols; k++) {
      int j = 0;
      for (int member : this.getColumn(k)) {
        if (j > 0) {
          colsWriter.write(" ");
        }
        j++;
        colsWriter.write(Integer.toString(member + 1)); // Add 1 to get back 1-indexing
      }
      colsWriter.newLine();
    }
    colsWriter.close();
  }

  /**
   * Get a row as a list of nonzero indices.
   *
   * @param rowId (required) row id
   * @return sparse list representation of row
   */
  public int[] getRow(int rowId) {
    return this.rows[rowId];
  }

  /**
   * Get a column as a set of nonzero indices.
   *
   * @param colId (required) column id
   * @return sparse set representation of column
   */
  public IntSet getColumn(int colId) {
    return this.cols[colId];
  }

  int[] getColumnAsArray(int colId) {
    IntSet col = this.getColumn(colId);
    int[] ret = new int[col.size()];
    IntIterator iter = col.iterator();
    for (int i = 0; iter.hasNext(); i++) {
      ret[i] = iter.nextInt();
    }
    return ret;
  }

  /**
   * Average number of nonzeros per row.
   *
   * @return average number of nonzeros per row
   */
  public double nnzPerRow() {
    double count = 0;
    for (int i = 0; i < this.numRows; i++) {
      count += this.rows[i].length;
    }
    return count / this.numRows;
  }

  /**
   * Ratio of empty rows among all rows.
   * Empty row means no assignment has been made.
   *
   * @return proportion of empty rows
   */
  public double emptyRowProportion() {
    double count = 0;
    for (int i = 0; i < this.numRows; i++) {
      if (this.rows[i].length == 0) {
        count++;
      }
    }
    return count / this.numRows;
  }

  double emptyColProportion() {
    double count = 0;
    for (int i = 0; i < this.numCols; i++) {
      if (this.cols[i].isEmpty()) {
        count++;
      }
    }
    return count / this.numCols;
  }

  /**
   * Size of column with fewest nonzeros, except for empty columns.
   *
   * @return smallest column size
   */
  int minColSizeAboveZero() {
    int min = Integer.MAX_VALUE;
    for (IntSet col : this.cols) {
      if (col.size() > 0 && col.size() < min) {
        min = col.size();
      }
    }
    return min;
  }

  /**
   * Size of column with the most nonzeros.
   *
   * @return biggest column size
   */
  int maxColSize() {
    int max = Integer.MIN_VALUE;
    for (IntSet col : this.cols) {
      if (col.size() > max) {
        max = col.size();
      }
    }
    return max;
  }


  /**
   * Set a row to a new value. Updates both rows and cols view.
   * Since rows are sorted, find intersection by 2-way merge.
   *
   * @param rowId (required) row id
   * @param newRow (required) value of the new row
   */
  void updateRow(int rowId, int[] newRow) {
    int oldHead = 0;
    int newHead = 0;
    int[] oldRow = this.rows[rowId];
    while (oldHead < oldRow.length && newHead < newRow.length) {
      int oldColId = oldRow[oldHead];
      int newColId = newRow[newHead];
      if (oldColId < newColId) {
        synchronized (this.getColumn(oldColId)) {
          this.getColumn(oldColId).remove(rowId);
        }
        oldHead++;
      } else if (oldColId > newColId) {
        synchronized (this.getColumn(newColId)) {
          this.getColumn(newColId).add(rowId);
        }
        newHead++;
      } else {
        oldHead++;
        newHead++;
      }
    }
    while (oldHead < oldRow.length) {
      int oldColId = oldRow[oldHead];
      synchronized (this.getColumn(oldColId)) {
        this.getColumn(oldColId).remove(rowId);
      }
      oldHead++;
    }
    while (newHead < newRow.length) {
      int newColId = newRow[newHead];
      synchronized (this.getColumn(newColId)) {
        this.getColumn(newColId).add(rowId);
      }
      newHead++;
    }
    this.rows[rowId] = newRow;
  }

  /**
   * Compute the number of nonzeros in the sum of given columns.
   * Equivalently, take union of the nonzero indices of columns and return its size.
   *
   * @param colIds (required) list of column id's to sum over
   * @return number of nonzeros in the column sum
   */
  int nnzInColumnSum(int[] colIds) {
    // Give a guess to the size of the union
    int totalSize = 0;
    for (int colId : colIds) {
      totalSize += this.getColumn(colId).size();
    }
    // Take union
    IntSet colIdUnion = new IntOpenHashSet(totalSize, 0.99F);
    for (int colId : colIds) {
      synchronized (this.getColumn(colId)) {
        colIdUnion.addAll(this.getColumn(colId));
      }
    }
    return colIdUnion.size();
  }

  /**
   * Neighborhood initialization using conductance.
   *
   * @param numNeighborhoods (required) initial number of communities
   * @return a list of sets, where each set represent an initial community
   */
  private IntSet[] getNeighborhoodInitialization(
      Graph g,
      int numNeighborhoods,
      BiFunction<Graph, Integer, Double> scoreForNeighborhoodFn) {
    System.err.println("Going to get scores for all vertices");
    double tic = System.currentTimeMillis();
    double[] conductanceScores = new double[this.numRows];
    Integer[] vertexIds = new Integer[this.numRows];
    for (int i = 0; i < this.numRows; i++) {
      conductanceScores[i] = scoreForNeighborhoodFn.apply(g, i);
      vertexIds[i] = i;
      if (i > 0 && i % 1000000 == 0) {
        System.err.format("Done getting scores for %d vertices\r", i);
      }
    }
    // Sort conductance scores in ascending order
    Arrays.sort(vertexIds,
        (i1, i2) -> (int) Math.signum(conductanceScores[i1] - conductanceScores[i2]));
    double toc = System.currentTimeMillis();
    System.err.println("Got scores for all vertices, time: " + (toc - tic) / 1000 + " secs");

    // Loop over from lowest conductance, assign each neighborhood to a community
    IntSet[] initCommunity = new IntSet[numNeighborhoods];
    IntSet seen = new IntOpenHashSet(this.numRows);
    int kSoFar = 0;
    for (int i = 0; i < this.numRows; i++) {
      int vertexId = vertexIds[i];
      if (seen.contains(vertexId)) {
        continue;
      }
      IntSet c = new IntOpenHashSet(g.getDegree(vertexId));
      for (int l = 0; l < g.getDegree(vertexId); l++) {
        c.add(g.getNeighbors(vertexId)[l]);
      }
      c.add(vertexId);
      seen.addAll(c);
      initCommunity[kSoFar++] = c;
      if (kSoFar >= numNeighborhoods) {
        break;
      }
    }
    // Assign random vertex to empty community
    while (kSoFar < numNeighborhoods) {
      IntSet randomVertexIds = new IntOpenHashSet();
      randomVertexIds.add(ThreadLocalRandom.current().nextInt(this.numRows));
      initCommunity[kSoFar++] = randomVertexIds;
    }
    return initCommunity;
  }

}
