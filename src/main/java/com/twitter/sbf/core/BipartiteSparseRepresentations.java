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

/**
 * Class that holds latent representation vectors for vertices of a bipartite graph with
 * community structure.
 */
public class BipartiteSparseRepresentations {
  public SparseRealMatrix leftRepresentations;
  public SparseRealMatrix rightRepresentations;

  /**
   * Constructor
   * @param leftRepresentations SparseRealMatrix object containing representation vectors for left
   * @param rightRepresentations SparseRealMatrix object containing representation vectors for right
   */
  public BipartiteSparseRepresentations(SparseRealMatrix leftRepresentations,
                                        SparseRealMatrix rightRepresentations) {
    if (leftRepresentations.getNumCols() != rightRepresentations.getNumCols()) {
      throw new RuntimeException(
          "Left and right representations must have the same dimensionality."
      );
    }
    this.leftRepresentations = leftRepresentations;
    this.rightRepresentations = rightRepresentations;
  }
}
