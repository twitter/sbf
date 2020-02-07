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

public class WeightedId {
  public int neighborId;
  public float weight;

  public WeightedId(int nId, float wt) {
    neighborId = nId;
    weight = wt;
  }

  public boolean equalsWithTolerance(Object o2) {
    WeightedId w = o2 instanceof WeightedId ? ((WeightedId) o2) : null;
    return w.neighborId == this.neighborId && Math.abs(w.weight - this.weight) < 1e-5;
  }
}
