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

package com.twitter.sbf.generator;

import com.twitter.sbf.graph.Graph;

import it.unimi.dsi.fastutil.ints.IntSet;

public class GraphAndGroundTruth {
  private Graph g;
  private IntSet[] clusters;

  GraphAndGroundTruth(Graph g, IntSet[] clusters) {
    this.g = g;
    this.clusters = clusters;
  }

  public Graph getGraph() {
    return g;
  }

  public IntSet[] getGroundTruth() {
    return clusters;
  }
}
