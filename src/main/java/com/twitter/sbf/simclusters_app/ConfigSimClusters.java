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

package com.twitter.sbf.simclusters_app;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomAdaptor;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

import com.twitter.sbf.core.AlgorithmConfig;
import com.twitter.sbf.core.ProposalStrategy;

/**
 * Stores all the configurations required by the SimClusters algorithm.
 */
public class ConfigSimClusters {
  private String metisFile = "metis_data"; // Bipartite graph stored in metis format
  //Number of communities/latent features
  private int k;
  //thresholding parameters
  private double thresholdForStepOne = 0.0;
  private int thresholdForStepThree = 1;
  //Whether top-k retention on edges should be done, and for what value of k
  private boolean retainTopK = false;
  private int topKParameter = 1;
  //Whether weights should be squared in projected Graph
  private boolean squareWeights = false;
  //Whether columns of right representations matrix should be normalized to unit norm.
  private boolean normalizeColumnsRight = false;
  //parameters for the second step of the algorithm
  private AlgorithmConfig mhAlgorithmParams = new AlgorithmConfig();
  private boolean initFromRandomNeighborhoods =  true;
  private boolean initFromBestNeighborhood = false;
  private boolean initFromNonoverlappingNeighborhood = false;
  //Whether rows of left representations matrix should be normalized to unit norm.
  private boolean normalizeRowsLeft = false;
  //Whether fourth step should be performed
  private boolean applyLogTransform = false;
  private boolean applyStepFour = false;
  //If step four is applied, what is the threshold
  private double thresholdForStepFour = 0.0;
  //Name of files where representations for left and right vertices need to printed.
  private String outputLeftRepFile = "left_rep";
  private String outputRightRepFile = "right_rep";


  /**
   * Constructor that reads from values from file
   * @param filename
   * @throws IOException
   */
  ConfigSimClusters(String filename) throws IOException {
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
          this.mhAlgorithmParams.cpu = Integer.parseInt(value);
          break;
        case "k":
          this.k = Integer.parseInt(value);
          this.mhAlgorithmParams.k = this.k;
          break;
        case "scaleCoeff":
          this.mhAlgorithmParams.scaleCoeff = Double.parseDouble(value);
          break;
        case "maxEpoch":
          this.mhAlgorithmParams.maxEpoch = Integer.parseInt(value);
          break;
        case "eps":
          this.mhAlgorithmParams.eps = Double.parseDouble(value);
          break;
        case "evalRatio":
          this.mhAlgorithmParams.evalRatio = Double.parseDouble(value);
          break;
        case "evalEvery":
          this.mhAlgorithmParams.evalEvery = Integer.parseInt(value);
          break;
        case "metisFile":
          this.metisFile = value;
          break;
        case "outputLeftRepFile":
          this.outputLeftRepFile = value;
          break;
        case "outputRightRepFile":
          this.outputRightRepFile = value;
          break;
        case "thresholdForStepOne":
          this.thresholdForStepOne = Double.parseDouble(value);
          break;
        case "thresholdForStepThree":
          this.thresholdForStepThree = Integer.parseInt(value);
          break;
        case "updateImmediately":
          this.mhAlgorithmParams.updateImmediately = Boolean.parseBoolean(value);
          break;
        case "noLocking":
          this.mhAlgorithmParams.noLocking = Boolean.parseBoolean(value);
          break;
        case "minClusterSize":
          this.mhAlgorithmParams.minClusterSize = Integer.parseInt(value);
          break;
        case "proposalStrategy":
          if ("FractionOfNeighbors".equalsIgnoreCase(value)) {
            this.mhAlgorithmParams.proposalStrategy
                = ProposalStrategy.FractionOfNeighbors;
          } else if ("FractionOfNeighborsAndCluster".equalsIgnoreCase(value)) {
            this.mhAlgorithmParams.proposalStrategy
                = ProposalStrategy.FractionOfNeighborsAndCluster;
          } else if ("SingleOrZeroMembershipLikelihood".equalsIgnoreCase(value)) {
            this.mhAlgorithmParams.proposalStrategy
                = ProposalStrategy.SingleOrZeroMembershipLikelihood;
          } else if ("MultipleMembershipLikelihood".equalsIgnoreCase(value)) {
            this.mhAlgorithmParams.proposalStrategy
                = ProposalStrategy.MultipleMembershipLikelihood;
          } else {
            System.err.println("Unrecognized value for proposalStrategy: "
                + value + ". Setting to default SingleOrZeroMembershipLikelihood");
            this.mhAlgorithmParams.proposalStrategy
                = ProposalStrategy.SingleOrZeroMembershipLikelihood;
          }
          break;
        case "maxMembershipsPerVertex":
          this.mhAlgorithmParams.maxMembershipsPerVertex = Integer.parseInt(value);
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
          this.mhAlgorithmParams.runAllEpochs = Boolean.parseBoolean(value);
          break;
        case "useWeightCoeffForProposal":
          this.mhAlgorithmParams.useWeightCoeffForProposal = Boolean.parseBoolean(value);
          break;
        case "randomSeed":
          int seed = Integer.parseInt(value);
          System.out.println("randomSeed = " + seed);
          this.mhAlgorithmParams.rng =
              new RandomAdaptor(new SynchronizedRandomGenerator(new JDKRandomGenerator(seed)));
          break;
        case "divideResultIntoConnectedComponents":
          this.mhAlgorithmParams.divideResultIntoConnectedComponents = Boolean.parseBoolean(value);
          break;
        case "removeWeakLinksForConnectedComponents":
          this.mhAlgorithmParams.removeWeakLinksForConnectedComponents
            = Boolean.parseBoolean(value);
          break;
        case "wtCoeff":
          this.mhAlgorithmParams.wtCoeff = Double.parseDouble(value);
          break;
        case "useTemperatureSchedule":
          this.mhAlgorithmParams.useTemperatureSchedule = Boolean.parseBoolean(value);
          break;
        case "maxTemperature":
          this.mhAlgorithmParams.maxTemperature = Double.parseDouble(value);
          break;
        case "retainTopK":
          this.retainTopK = Boolean.parseBoolean(value);
          break;
        case "topKParameter":
          this.topKParameter = Integer.parseInt(value);
          break;
        case "squareWeights":
          this.squareWeights = Boolean.parseBoolean(value);
          break;
        case "normalizeColumnsRight":
          this.normalizeColumnsRight = Boolean.parseBoolean(value);
          break;
        case "normalizeRowsLeft":
          this.normalizeRowsLeft = Boolean.parseBoolean(value);
          break;
        case "applyStepFour":
          this.applyStepFour = Boolean.parseBoolean(value);
          break;
        case "thresholdForStepFour":
          this.thresholdForStepFour = Double.parseDouble(value);
          break;
        case "applyLogTransform":
          this.applyLogTransform = Boolean.parseBoolean(value);
          break;
        default:
          System.err.println("Unknown config name: " + key);
          break;
      }
    }
    br.close();
  }
  //Accessor methods
  String getMetisFile() {
    return this.metisFile;
  }

  AlgorithmConfig getAlgoConfig() {
    return this.mhAlgorithmParams;
  }

  int getK() {
    return this.k;
  }

  double getThresholdForStepOne() {
    return this.thresholdForStepOne;
  }

  int getThresholdForStepThree() {
    return this.thresholdForStepThree;
  }

  boolean getSquareWeights() {
    return this.squareWeights;
  }

  boolean getNormalizeColumnsRight() {
    return this.normalizeColumnsRight;
  }

  boolean getRetainTopK() {
    return this.retainTopK;
  }

  int getTopKParameter() {
    return this.topKParameter;
  }

  boolean getInitFromRandomNeighborhoods() {
    return this.initFromRandomNeighborhoods;
  }

  boolean getInitFromBestNeighborhood() {
    return this.initFromBestNeighborhood;
  }

  boolean getInitFromNonoverlappingNeighborhood() {
    return this.initFromNonoverlappingNeighborhood;
  }

  boolean getNormalizeRowsLeft() {
    return this.normalizeRowsLeft;
  }

  boolean getApplyLogTransform() {
    return this.applyLogTransform;
  }

  boolean getApplyStepFour() {
    return this.applyStepFour;
  }

  double getThresholdForStepFour() {
    return this.thresholdForStepFour;
  }

  String getOutputLeftRepFile() {
    return this.outputLeftRepFile;
  }

  String getOutputRightRepFile() {
    return this.outputRightRepFile;
  }
}
