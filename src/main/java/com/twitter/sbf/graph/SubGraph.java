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

import java.util.Set;

import com.twitter.sbf.util.IntIteratorFromArray;
import com.twitter.sbf.util.Util;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntList;

public class SubGraph implements GraphInterface {
  private Graph fullGraph;
  private Set<Integer> nodesOfSubgraph;
  public SubGraph(Graph fullGraph, Set<Integer> nodesOfSubgraph) {
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
      IntList neighbors = new IntArrayList();
      for (int n : allNeighbors) {
        if (nodesOfSubgraph.contains(n)) {
          neighbors.add(n);
        }
      }
      return new IntIteratorFromArray(neighbors.toIntArray());
    }
  }
}
