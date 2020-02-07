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

package com.twitter.sbf.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

/**
 * Utility class for the package.
 * <p>
 */
public final class Util {
  private Util() {
  }

  /**
   * Randomly choose k numbers from [0, 1, ..., n-1]
   *
   * @param n (required) Size of the pool
   * @param k (required) Size of the selection
   * @return An ArrayList of k numbers
   */
  public static int[] reservoirSampling(int n, int k) {
    int[] result = new int[k];
    if (k == 0) {
      return result;
    }
    for (int i = 0; i < k; i++) {
      result[i] = i;
    }
    for (int i = k; i < n; i++) {
      int j = ThreadLocalRandom.current().nextInt(i + 1);
      if (j < k) {
        result[j] = i;
      }
    }
    return result;
  }

  private static void testReservoirSampling() {
    int n = 10;
    int k = 3;
    int repeat = 10000;
    int[] freq = new int[n];
    Arrays.fill(freq, 0);
    for (int i = 0; i < repeat; i++) {
      int[] result = Util.reservoirSampling(n, k);
      for (int r : result) {
        freq[r] += 1;
      }
    }
    for (int i = 0; i < n; i++) {
      System.out.println("index: " + i + " prob: " + ((double) (freq[i]) / repeat));
    }
  }

  /**
   * Test if two sorted lists have at least one common element.
   * Does NOT check if input is sorted.
   *
   * @param a1 (required) a sorted list
   * @param a2 (required) a sorted list
   * @return true iff two sorted lists have at least one element in common
   */
  public static boolean hasCommonElement(final int[] a1, final int[] a2) {
    int head1 = 0;
    int head2 = 0;
    while (head1 < a1.length && head2 < a2.length) {
      int value1 = a1[head1];
      int value2 = a2[head2];
      if (value1 < value2) {
        head1++;
      } else if (value1 > value2) {
        head2++;
      } else {
        return true;
      }
    }
    return false;
  }

  /**
   * Split an array into some chunks and return the ranges.
   *
   * @param arraySize size of the array
   * @param maxNumChunks max number of chunks to split the array into.
   * @return array of ranges of each chunk
   */
  public static Range[] chunkArray(int arraySize, int maxNumChunks) {
    int maxChunkSize = (arraySize % maxNumChunks == 0)
        ? arraySize / maxNumChunks
        : arraySize / maxNumChunks + 1;
    int numChunks;
    if (maxChunkSize == 0) {
      numChunks = 1;
    } else {
      numChunks = (arraySize % maxChunkSize == 0)
          ? arraySize / maxChunkSize
          : arraySize / maxChunkSize + 1;
    }

    Range[] ranges = new Range[numChunks];
    for (int i = 0; i < numChunks; i++) {
      int start = maxChunkSize * i;
      int end = Math.min(maxChunkSize * (i + 1), arraySize);
      ranges[i] = new Range(start, end);
    }
    return ranges;
  }

  /**
   * Repeat a string n times.
   *
   * @param n (required) num of repeat
   * @param s (required) string to repeat
   * @return a new string with s repeating n times
   */
  public static String repeatString(int n, String s) {
    return new String(new char[n]).replace("\0", s);
  }

  /**
   * Count the number of lines of a file.
   *
   * @param filename (required) file path
   * @return number of lines
   */
  public static long getNumLinesOfFile(String filename) throws IOException {
    return Files.lines(Paths.get(filename)).count();
  }

  /**
   * Check if a file already exists.
   *
   * @param filename (required) file to test existance
   * @return true iff the file exists
   */
  public static boolean fileExists(String filename) {
    File f = new File(filename);
    return f.isFile();
  }

  /**
   * Append values as a space-separated line in a file.
   *
   * @param filename (required) filename to append to
   * @param values (required) values to append
   */
  public static void appendDoublesToFile(String filename, double[] values) throws IOException {
    BufferedWriter bw = new BufferedWriter(new FileWriter(filename, true));
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        bw.write("\n");
      }
      bw.write(String.format("%.2f", values[i]));
    }
    bw.newLine();
    bw.flush();
    bw.close();
  }

  /**
   * Convert set of integers to array
   *
   * @param set set to be converted
   * @return array
   */
  public static int[] intArrayFromSet(Set<Integer> set) {
    int[] ret = new int[set.size()];
    int index = 0;
    for (Integer i : set) {
      ret[index++] = i;
    }
    return ret;
  }

  /**
   * convert iterator to set
   * @param iter input iterator
   * @return set
   */
  public static IntSet setFromIntIterator(IntIterator iter) {
    IntSet ret = new IntOpenHashSet();
    while (iter.hasNext()) {
      ret.add(iter.nextInt());
    }
    return ret;
  }
}
