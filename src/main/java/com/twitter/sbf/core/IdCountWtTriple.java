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

public class IdCountWtTriple {
  public int id;
  public int count;
  public double wt;

  IdCountWtTriple(int id, int count, double wt) {
    this.id = id;
    this.count = count;
    this.wt = wt;
  }

  @Override
  public String toString() {
    return String.format("id %d, count %d, wt %g", id, count, wt);
  }
}
