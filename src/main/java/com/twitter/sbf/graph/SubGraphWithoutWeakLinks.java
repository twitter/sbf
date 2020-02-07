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

import java.util.Arrays;
import java.util.Set;

import com.twitter.sbf.util.IntIteratorFromArray;
import com.twitter.sbf.util.Util;

import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

/**
 * An edge (u, v) is called a "weak link" if u and v share no common neighbors.
 *
 * This class creates a representation of a subgraph of an input graph with 2 properties:
 * 1. The caller can specify which nodes of the original graph are included in the subgraph.
 * 2. Only edges which are not weak links are included in the neighbors returned for a node.
 */
public class SubGraphWithoutWeakLinks implements GraphInterface {
  private Graph fullGraph;
  private Set<Integer> nodesOfSubgraph;
  private long numWeakLinks = 0;

  public SubGraphWithoutWeakLinks(Graph fullGraph, Set<Integer> nodesOfSubgraph) {
    this.fullGraph = fullGraph;
    this.nodesOfSubgraph = nodesOfSubgraph;
  }

  /**
   * Return the ids of the nodes in this graph
   * @return iterator with the ids
   */
  public IntIterator getAllNodeIds() {
    return new IntIteratorFromArray(Util.intArrayFromSet(nodesOfSubgraph));
  }

  /**
   * Return the neighbors of a specified node
   * @param nodeId id of the node whose neighbors are desired
   * @return iterator with the neighbors of the specified node
   */
  public IntIterator getNeighbors(int nodeId) {
    if (!nodesOfSubgraph.contains(nodeId)) {
      return IntIterators.EMPTY_ITERATOR;
    } else {
      int[] allNeighbors = fullGraph.getNeighbors(nodeId);
      IntSet neighborsIncludingWeakLinks = new IntOpenHashSet();
      for (int n : allNeighbors) {
        if (nodesOfSubgraph.contains(n)) {
          neighborsIncludingWeakLinks.add(n);
        }
      }

      int[] neighborsWithoutWeakLinks = new int[neighborsIncludingWeakLinks.size()];
      int lengthOfNeighborsWithoutWeakLinks = 0;
      for (int neighbor : neighborsIncludingWeakLinks) {
        boolean sharedNeighborsSoFar = false;
        for (int neighborOfNeighbor: fullGraph.getNeighbors(neighbor)) {
          if (neighborsIncludingWeakLinks.contains(neighborOfNeighbor)
              && nodesOfSubgraph.contains(neighborOfNeighbor)) {
            sharedNeighborsSoFar = true;
            break;
          }
        }
        if (sharedNeighborsSoFar) {
          neighborsWithoutWeakLinks[lengthOfNeighborsWithoutWeakLinks++] = neighbor;
        } else {
          numWeakLinks++;
        }

      }

      int[] ret = Arrays.copyOf(neighborsWithoutWeakLinks, lengthOfNeighborsWithoutWeakLinks);

      return new IntIteratorFromArray(ret);
    }
  }

  public long getNumWeakLinks() {
    return numWeakLinks;
  }
}
