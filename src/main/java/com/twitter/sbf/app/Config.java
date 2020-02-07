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

package com.twitter.sbf.app;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import com.google.common.collect.ImmutableList;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomAdaptor;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

import com.twitter.sbf.core.AlgorithmConfig;
import com.twitter.sbf.core.ProposalStrategy;

/**
 * Stores all the configurations required by the algorithm.
 */
public class Config {
  String metisFile = ""; // graph stored in metis format
  private AlgorithmConfig algoConfig = new AlgorithmConfig();

  boolean initFromRandomNeighborhoods = true;
  boolean initFromBestNeighborhood = false; // initialize from neighborhoods with best conductance
  boolean initFromNonoverlappingNeighborhood = false;
  String initFromRowsFile = ""; // file where each line is a row of Z
  String initFromColsFile = ""; // file where each line is a column of Z
  String outputByRowsFile = ""; // one line per vertex, indicating all clusters for this vertex
  String outputByColsFile = ""; // one line per cluster, indicating all members of that community
  String clusterPrecisionFile = ""; // file where every line is a cluster precision
  String outputRowsWithScoresFile = ""; // one line per vertex, giving clusters along with scores

  Config(String filename) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(filename));
    String line;
    while ((line = br.readLine()) != null) {
      // Ignore empty line or comments starting with %
      if (line.isEmpty() || line.startsWith("%")) {
        continue;
      }
      String[] tokens = line.trim().split("\\s+");
      assert tokens.length == 2;
      String key = tokens[0];
      String value = tokens[1];
      switch (key) {
        case "cpu":
          this.algoConfig.cpu = Integer.parseInt(value);
          break;
        case "K":
          this.algoConfig.k = Integer.parseInt(value);
          break;
        case "scaleCoeff":
          this.algoConfig.scaleCoeff = Double.parseDouble(value);
          break;
        case "maxEpoch":
          this.algoConfig.maxEpoch = Integer.parseInt(value);
          break;
        case "eps":
          this.algoConfig.eps = Double.parseDouble(value);
          break;
        case "evalRatio":
          this.algoConfig.evalRatio = Double.parseDouble(value);
          break;
        case "evalEvery":
          this.algoConfig.evalEvery = Integer.parseInt(value);
          break;
        case "metisFile":
          this.metisFile = value;
          break;
        case "initFromRowsFile":
          this.initFromRowsFile = value;
          break;
        case "initFromColsFile":
          this.initFromColsFile = value;
          break;
        case "outputByRowsFile":
          this.outputByRowsFile = value;
          break;
        case "outputByColsFile":
          this.outputByColsFile = value;
          break;
        case "outputRowsWithScoresFile":
          this.outputRowsWithScoresFile = value;
          break;
        case "clusterPrecisionFile":
          this.clusterPrecisionFile = value;
          break;
        case "updateImmediately":
          this.algoConfig.updateImmediately = Boolean.parseBoolean(value);
          break;
        case "noLocking":
          this.algoConfig.noLocking = Boolean.parseBoolean(value);
          break;
        case "minClusterSize":
          this.algoConfig.minClusterSize = Integer.parseInt(value);
          break;
        case "proposalStrategy":
          if ("FractionOfNeighbors".equalsIgnoreCase(value)) {
            this.algoConfig.proposalStrategy = ProposalStrategy.FractionOfNeighbors;
          } else if ("FractionOfNeighborsAndCluster".equalsIgnoreCase(value)) {
            this.algoConfig.proposalStrategy = ProposalStrategy.FractionOfNeighborsAndCluster;
          } else if ("SingleOrZeroMembershipLikelihood".equalsIgnoreCase(value)) {
            this.algoConfig.proposalStrategy = ProposalStrategy.SingleOrZeroMembershipLikelihood;
          } else if ("MultipleMembershipLikelihood".equalsIgnoreCase(value)) {
            this.algoConfig.proposalStrategy = ProposalStrategy.MultipleMembershipLikelihood;
          } else {
            System.err.println("Unrecognized value for proposalStrategy: "
                + value + ". Setting to default SingleOrZeroMembershipLikelihood");
            this.algoConfig.proposalStrategy = ProposalStrategy.SingleOrZeroMembershipLikelihood;
          }
          break;
        case "maxMembershipsPerVertex":
          this.algoConfig.maxMembershipsPerVertex = Integer.parseInt(value);
          break;
        case "initFromNonoverlappingNeighborhood":
          this.initFromNonoverlappingNeighborhood = Boolean.parseBoolean(value);
          break;
        case "initFromRandomNeighborhood":
          this.initFromRandomNeighborhoods = Boolean.parseBoolean(value);
          break;
        case "initFromBestNeighborhood":
          this.initFromBestNeighborhood = Boolean.parseBoolean(value);
          break;
        case "runAllEpochs":
          this.algoConfig.runAllEpochs = Boolean.parseBoolean(value);
          break;
        case "useWeightCoeffForProposal":
          this.algoConfig.useWeightCoeffForProposal = Boolean.parseBoolean(value);
          break;
        case "randomSeed":
          int seed = Integer.parseInt(value);
          System.out.println("randomSeed = " + seed);
          this.algoConfig.rng =
              new RandomAdaptor(new SynchronizedRandomGenerator(new JDKRandomGenerator(seed)));
          break;
        case "divideResultIntoConnectedComponents":
          this.algoConfig.divideResultIntoConnectedComponents = Boolean.parseBoolean(value);
          break;
        case "removeWeakLinksForConnectedComponents":
          this.algoConfig.removeWeakLinksForConnectedComponents = Boolean.parseBoolean(value);
          break;
        case "wtCoeff":
          this.algoConfig.wtCoeff = Double.parseDouble(value);
          break;
        case "useTemperatureSchedule":
          this.algoConfig.useTemperatureSchedule = Boolean.parseBoolean(value);
          break;
        case "maxTemperature":
          this.algoConfig.maxTemperature = Double.parseDouble(value);
          break;
        default:
          System.err.println("Unknown config name: " + key);
          break;
      }
    }
    br.close();
  }

  @Override
  public String toString() {
    return String.join(
        "\n",
        ImmutableList.of(
            "metisFile = " + metisFile,
            "initFromRowsFile = " + initFromRowsFile,
            "initFromColsFile = " + initFromColsFile,
            "initFromRandomNeighborhoods = " + initFromRandomNeighborhoods,
            "initFromBestNeighborhood = " + initFromBestNeighborhood,
            "outputByRowsFile = " + outputByRowsFile,
            "outputByColsFile = " + outputByColsFile,
            "clusterPrecisionFile = " + clusterPrecisionFile,
            "outputRowsWithScoresFile = " + outputRowsWithScoresFile,
            algoConfig.toString()
        )
    );
  }


  AlgorithmConfig getAlgoConfig() {
    return this.algoConfig;
  }
}
