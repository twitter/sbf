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

import com.google.common.collect.ImmutableList;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomAdaptor;

public class AlgorithmConfig {
  public int cpu = Runtime.getRuntime().availableProcessors(); // size of the thread pool
  public int k = 100; // number of communities
  public double scaleCoeff = 0.5; // how much log-likelihood is scaled
  public int maxEpoch = 500; // max number of sweep over the vertices
  // do not terminate based on convergence criteria; run for all epochs
  public boolean runAllEpochs = false;
  public double eps = 1e-2; // stop if accept rate is smaller than this
  public double evalRatio = 0.001; // proportion of nodesOfSubgraph to evaluate metrics, saves time
  public int evalEvery = 1; // evaluate metrics every X epochs
  public boolean useWeightCoeffForProposal = false; // use weight coeff in proposal distribution
  public double maxTemperature = 100;
  public double temperatureRatio = 0.95;
  public boolean useTemperatureSchedule = false;
  public RandomAdaptor rng = new RandomAdaptor(new JDKRandomGenerator());
  public double wtCoeff = 5.0;
  public boolean divideResultIntoConnectedComponents = true;
  public boolean removeWeakLinksForConnectedComponents = true;
  public boolean updateImmediately = true;
  public boolean noLocking = true;
  public ProposalStrategy proposalStrategy = ProposalStrategy.SingleOrZeroMembershipLikelihood;
  public int maxMembershipsPerVertex = 5;
  public int minClusterSize = 2;
  public boolean ignoreSameProposal = true;

  @Override
  public String toString() {
    return String.join(
        "\n",
        ImmutableList.of(
            "cpu = " + cpu,
            "K = " + k,
            "scaleCoeff = " + scaleCoeff,
            "maxEpoch = " + maxEpoch,
            "runAllEpochs = " + runAllEpochs,
            "eps = " + eps,
            "evalRatio = " + evalRatio,
            "evalEvery = " + evalEvery,
            "useWeightCoeffForProposal = " + useWeightCoeffForProposal,
            "maxTemperature = " + maxTemperature,
            "temperatureRatio = " + temperatureRatio,
            "wtCoeff = " + wtCoeff,
            "divideResultIntoConnectedComponents = " + divideResultIntoConnectedComponents,
            "removeWeakLinksForConnectedComponents = " + removeWeakLinksForConnectedComponents,
            "updateImmediately = " + updateImmediately,
            "noLocking = " + noLocking,
            "proposalStrategy = " + proposalStrategy,
            "maxMembershipsPerVertex = " + maxMembershipsPerVertex,
            "minClusterSize = " + minClusterSize
        )
    );
  }

  public AlgorithmConfig withK(int k2) {
    this.k = k2;
    return this;
  }

  public AlgorithmConfig withMaxEpoch(int me) {
    this.maxEpoch = me;
    return this;
  }

  public AlgorithmConfig withEvalRatio(double er) {
    this.evalRatio = er;
    return this;
  }

  public AlgorithmConfig withRunAllEpochs(boolean runAll) {
    this.runAllEpochs = runAll;
    return this;
  }

  public AlgorithmConfig withCpu(int c) {
    this.cpu = c;
    return this;
  }

  public AlgorithmConfig withScaleCoeff(double sc) {
    this.scaleCoeff = sc;
    return this;
  }

  public AlgorithmConfig withUseWeightCoeffForProposal(boolean b) {
    this.useWeightCoeffForProposal = b;
    return this;
  }

  public AlgorithmConfig withMaxTemperature(int temp) {
    this.maxTemperature = temp;
    return this;
  }

  public AlgorithmConfig withUseTemperatureSchedule(boolean b) {
    this.useTemperatureSchedule = b;
    return this;
  }

  public AlgorithmConfig withRng(RandomAdaptor r) {
    this.rng = r;
    return this;
  }

  public AlgorithmConfig withWtCoeff(double dm) {
    this.wtCoeff = dm;
    return this;
  }
}
