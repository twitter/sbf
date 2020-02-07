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

import com.twitter.sbf.util.IntIteratorFromArray;

import it.unimi.dsi.fastutil.ints.IntIterator;

public class GraphAdaptorForInterface implements GraphInterface {
  private Graph graph;
  public GraphAdaptorForInterface(Graph g) {
    graph = g;
  }

  public IntIterator getAllNodeIds() {
    return new IntIteratorFromArray(graph.getAllVertexIds());
  }

  public IntIterator getNeighbors(int nodeId) {
    return new IntIteratorFromArray(graph.getNeighbors(nodeId));
  }
}
