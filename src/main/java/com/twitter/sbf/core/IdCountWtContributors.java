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

import java.util.stream.Collectors;

import it.unimi.dsi.fastutil.ints.IntSet;

public class IdCountWtContributors {
  public int id;
  public int count;
  public double wt;
  public IntSet contributors;

  IdCountWtContributors(int id, int count, double wt, IntSet contributors) {
    this.id = id;
    this.count = count;
    this.wt = wt;
    this.contributors = contributors;
  }

  @Override
  public String toString() {
    return String.format("id %d, count %d, wt %g, contributors [%s]", id, count, wt,
        String.join(",",
            contributors.stream().map(i -> i.toString()).collect(Collectors.toList())));
  }
}
