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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import com.twitter.sbf.graph.ConnectedComponents;
import com.twitter.sbf.graph.Graph;
import com.twitter.sbf.util.SimpleIteratorFromIterator;

import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import static org.junit.Assert.assertTrue;

public class ConnectedComponentsTest {

  static Graph connected8NodeGraph() {
    List<String> lines = new ArrayList<>();
    lines.add("8 13");
    lines.add("2 3 4 5");
    lines.add("1 3 4");
    lines.add("1 2 4");
    lines.add("1 2 3");
    lines.add("1 6 7 8");
    lines.add("5 7 8");
    lines.add("5 6 8");
    lines.add("5 6 7");

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  private static Graph disConnected8NodeGraph() {
    List<String> lines = new ArrayList<>();
    lines.add("8 12");
    lines.add("2 3 4");
    lines.add("1 3 4");
    lines.add("1 2 4");
    lines.add("1 2 3");
    lines.add("6 7 8");
    lines.add("5 7 8");
    lines.add("5 6 8");
    lines.add("5 6 7");

    return new Graph(new SimpleIteratorFromIterator<>(lines.iterator()));
  }

  @Test
  public void test8NodeGraphs() {
    Set<IntSet> components = ConnectedComponents.connectedComponents(connected8NodeGraph());
    assertTrue(components.size() == 1);
    assertTrue(components.iterator().next().size() == 8);

    Set<IntSet> components2 = ConnectedComponents.connectedComponents(disConnected8NodeGraph());
    assertTrue(components2.size() == 2);
    Iterator<IntSet> iter = components2.iterator();
    assertTrue(iter.hasNext());
    Set<Integer> comp1 = iter.next();
    assertTrue(iter.hasNext());
    Set<Integer> comp2 = iter.next();

    assertTrue(comp1.size() == 4);
    assertTrue(comp2.size() == 4);
  }

  @Test
  public void testSubGraphConnectedComponents() {
    //Set<Integer> selectedNodes = new HashSet<>();
    IntSet selectedNodes = new IntOpenHashSet();
    selectedNodes.add(2);
    selectedNodes.add(3);
    selectedNodes.add(5);
    selectedNodes.add(6);


    Set<IntSet> comps =
        ConnectedComponents.subGraphConnectedComponents(connected8NodeGraph(), selectedNodes);
    assertTrue(comps.size() == 2);
    Iterator<IntSet> iter = comps.iterator();
    Set<Integer> comp1 = iter.next();
    Set<Integer> comp2 = iter.next();
    assertTrue(comp1.size() == 2 && comp2.size() == 2);
    assertTrue(
        comp1.contains(2) && comp1.contains(3) || (comp2.contains(2) && comp2.contains(3))
    );
    assertTrue(
        comp1.contains(5) && comp1.contains(6) || (comp2.contains(5) && comp2.contains(6))
    );
  }

  @Test
  public void testWeakLinkConnectedComponents() {
    Graph g = connected8NodeGraph();
    IntSet allNodes = new IntArraySet(g.getAllVertexIds());

    Set<IntSet> comps =
        ConnectedComponents.subGraphConnectedComponentsNoWeakLinks(g, allNodes);

    assertTrue(comps.size() == 2);
    Iterator<IntSet> iter = comps.iterator();
    Set<Integer> comp1 = iter.next();
    Set<Integer> comp2 = iter.next();
    Set<Integer> compWith0;
    Set<Integer> compWithout0;

    if (comp1.contains(0)) {
      compWith0 = comp1;
      compWithout0 = comp2;
    } else {
      compWith0 = comp2;
      compWithout0 = comp1;
    }

    assertTrue(compWith0.size() == 4);
    assertTrue(compWithout0.size() == 4);

    assertTrue(compWith0.contains(0) && compWith0.contains(1)
        && compWith0.contains(2) && compWith0.contains(3));
    assertTrue(compWithout0.contains(4) && compWithout0.contains(5)
        && compWithout0.contains(6) && compWithout0.contains(7));
  }

}
