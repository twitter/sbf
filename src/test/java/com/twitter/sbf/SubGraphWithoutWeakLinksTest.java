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

package com.twitter.sbf;

import org.junit.Test;

import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.graph.SubGraphWithoutWeakLinks;
import com.twitter.sbf.util.Util;

import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;

import static org.junit.Assert.assertTrue;

public class SubGraphWithoutWeakLinksTest {

  @Test
  public void testSubGraphWithoutWeaklinks() {
    Graph g = ConnectedComponentsTest.connected8NodeGraph();
    IntSet allNodes = new IntArraySet(g.getAllVertexIds());
    SubGraphWithoutWeakLinks sg = new SubGraphWithoutWeakLinks(g, allNodes);

    IntSet setFor0 = Util.setFromIntIterator(sg.getNeighbors(0));
    assertTrue(setFor0.size() == 3);
    assertTrue(setFor0.contains(1) && setFor0.contains(2) && setFor0.contains(3));

    IntSet setFor1 = Util.setFromIntIterator(sg.getNeighbors(1));
    assertTrue(setFor1.size() == 3);
    assertTrue(setFor1.contains(0) && setFor1.contains(2) && setFor1.contains(3));

    IntSet setFor2 = Util.setFromIntIterator(sg.getNeighbors(2));
    assertTrue(setFor2.size() == 3);
    assertTrue(setFor2.contains(0) && setFor2.contains(1) && setFor2.contains(3));

    IntSet setFor3 = Util.setFromIntIterator(sg.getNeighbors(3));
    assertTrue(setFor3.size() == 3);
    assertTrue(setFor3.contains(0) && setFor3.contains(1) && setFor3.contains(2));

    IntSet setFor4 = Util.setFromIntIterator(sg.getNeighbors(4));
    assertTrue(setFor4.size() == 3);
    assertTrue(setFor4.contains(5) && setFor4.contains(6) && setFor4.contains(7));

    IntSet setFor5 = Util.setFromIntIterator(sg.getNeighbors(5));
    assertTrue(setFor5.size() == 3);
    assertTrue(setFor5.contains(4) && setFor5.contains(6) && setFor5.contains(7));

    IntSet setFor6 = Util.setFromIntIterator(sg.getNeighbors(6));
    assertTrue(setFor6.size() == 3);
    assertTrue(setFor6.contains(5) && setFor6.contains(4) && setFor6.contains(7));

    IntSet setFor7 = Util.setFromIntIterator(sg.getNeighbors(7));
    assertTrue(setFor7.size() == 3);
    assertTrue(setFor7.contains(5) && setFor7.contains(6) && setFor7.contains(4));

  }
}
