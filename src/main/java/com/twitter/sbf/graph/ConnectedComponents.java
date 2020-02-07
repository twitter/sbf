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

import java.util.HashSet;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;

import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

/**
 * Utility methods for dividing a graph into its connected components.
 */
public final class ConnectedComponents {
  private ConnectedComponents() {
  }

  /**
   * Divide a graph into connected components
   * @param g Input graph
   * @return set of sets, each of which is a connected component
   */
  public static Set<IntSet> connectedComponents(Graph g) {
    return connectedComponents(new GraphAdaptorForInterface(g));
  }

  /**
   * Get the connected components of a graph
   * @param g input graph whose connected components are desired
   * @return set containing each connected component (each component is represented by an IntSet)
   */
  public static Set<IntSet> connectedComponents(GraphInterface g) {
    Set<IntSet> ret = new HashSet<>();
    Queue<Integer> q = new LinkedBlockingQueue<>();
    IntSet unSeen = new IntOpenHashSet();
    IntIterator nodeIter = g.getAllNodeIds();
    while (nodeIter.hasNext()) {
      unSeen.add(nodeIter.nextInt());
    }

    while (!unSeen.isEmpty()) {
      int newNode = unSeen.iterator().next();
      IntSet newComponent = new IntOpenHashSet();
      newComponent.add(newNode);
      q.add(newNode);
      unSeen.remove(newNode);
      while (!q.isEmpty()) {
        int currentNode = q.remove();
        IntIterator iter = g.getNeighbors(currentNode);
        while (iter.hasNext()) {
          int nbr = iter.nextInt();
          if (unSeen.contains(nbr)) {
            q.add(nbr);
            assert !newComponent.contains(nbr);
            newComponent.add(nbr);
            unSeen.remove(nbr);
          }
        }
      }
      ret.add(newComponent);
    }

    return ret;
  }

  /**
   * Connected components of a subgraph of a graph
   * @param g Input graph
   * @param nodesOfSubGraph nodes of subgraph
   * @return set of sets, each of which is a connected component of the subgraph
   */
  public static Set<IntSet> subGraphConnectedComponents(Graph g, Set<Integer> nodesOfSubGraph) {
    return connectedComponents(new SubGraph(g, nodesOfSubGraph));
  }

  /**
   * An edge (u, v) is called a "weak link" if u and v share no common neighbors.
   *
   * Calling this method will return the connected components of the subgraph of the input graph
   * after removing all the weak links from the subgraph.
   * @param g
   * @param nodesOfSubGraph
   * @return
   */
  public static Set<IntSet> subGraphConnectedComponentsNoWeakLinks(
      Graph g,
      Set<Integer> nodesOfSubGraph
  ) {
    return connectedComponents(new SubGraphWithoutWeakLinks(g, nodesOfSubGraph));
  }

}
