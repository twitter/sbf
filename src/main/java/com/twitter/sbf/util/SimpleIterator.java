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

import java.util.Optional;

/**
 * Simpler version of iterator. Had to create this interface instead of
 * reusing java.util.Iterator because it's hard to wrap BufferedReader.readLine() inside
 * java.util.Iterator (possible to implement using locks and what not, but it's annoying)
 */
public interface SimpleIterator<T> {
  /**
   * If this iterator has one more element, returns that element as the value of the Optional.
   * If not, it returns Optional.empty().
   * The caller needs to infer on the basis of the empty Optional that the iterator is exhausted.
   *
   * @return Either an element or empty.
   */
  Optional<T> next();
}
