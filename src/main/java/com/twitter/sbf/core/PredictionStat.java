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

public class PredictionStat {
  private double numActualPositive = 0;
  private double weightActualPositive = 0;
  private double numPredictedPositive = 0;
  private double weightPredictedPositive = 0;
  private double numTruePositive = 0;
  private double weightTruePositive = 0;
  private double numTrueNegative = 0;
  private double numEval = 0;
  private double numEvalVertices = 0;
  private double numEvalVerticesWithZeroTruePos = 0;

  void incActualPositive() {
    numActualPositive++;
  }

  void incPredictedPositive() {
    numPredictedPositive++;
  }

  void incTruePositive() {
    numTruePositive++;
  }

  void incTrueNegative() {
    numTrueNegative++;
  }

  void incEval(int add) {
    numEval += add;
  }

  void incEvalVertices() {
    numEvalVertices++;
  }

  void incEvalVerticesWithZeroTruePos() {
    numEvalVerticesWithZeroTruePos++;
  }

  void incWeightActualPositive(double add) {
    weightActualPositive += add;
  }

  void incWeightPredictedPositive(double add) {
    weightPredictedPositive += add;
  }

  void incWeightTruePositive(double add) {
    weightTruePositive += add;
  }

  /**
   * Add a stat to myself.
   *
   * @param other (required) the other stat object
   */
  public void add(PredictionStat other) {
    this.numActualPositive += other.numActualPositive;
    this.weightActualPositive += other.weightActualPositive;
    this.numPredictedPositive += other.numPredictedPositive;
    this.weightPredictedPositive += other.weightPredictedPositive;
    this.numTruePositive += other.numTruePositive;
    this.weightTruePositive += other.weightTruePositive;
    this.numTrueNegative += other.numTrueNegative;
    this.numEval += other.numEval;
    this.numEvalVertices += other.numEvalVertices;
    this.numEvalVerticesWithZeroTruePos += other.numEvalVerticesWithZeroTruePos;
  }

  public double precision() {
    return numTruePositive / numPredictedPositive;
  }

  double weightedPrecision() {
    if (Double.isNaN(weightPredictedPositive)) {
      return 0;
    } else {
      return weightTruePositive / weightPredictedPositive;
    }
  }

  double recall() {
    return numTruePositive / numActualPositive;
  }

  double weightedRecall() {
    if (weightActualPositive == 0) {
      return 0;
    } else {
      return weightTruePositive / weightActualPositive;
    }
  }

  /**
   * Harmonic mean of precision and recall
   */
  public double f1() {
    double prec = precision();
    double rec = recall();
    if (prec == 0 && rec == 0) {
      return 0;
    } else {
      return 2 * prec * rec / (prec + rec);
    }
  }

  double orphanRate() {
    return numEvalVerticesWithZeroTruePos / numEvalVertices;
  }

  /**
   * Harmonic mean of weighted precision and weighted recall
   */
  public double weightedF1() {
    double prec = weightedPrecision();
    double rec = weightedRecall();
    if (prec == 0 && rec == 0) {
      return 0;
    } else {
      return 2 * prec * rec / (prec + rec);
    }
  }

  /**
   * string representation, meant for logging
   */
  public String toString(boolean includeWeights) {
    String line1 = String.format(
        "Precision: %.2f, Recall: %.2f, F1: %.2f", precision(), recall(), f1()
    );
    String line2 = String.format(
        "Evals: %.0f, Actual Edges: %.0f, Predicted Edges: %.0f, True Positive: %.0f, "
        + "True Negative: %.0f", numEval, numActualPositive, numPredictedPositive, numTruePositive,
        numTrueNegative
    );
    String line3 = "";
    String line4 = "";
    if (includeWeights) {
      line3 = String.format(
          "\nActual Edges Weight: %.2f, Predicted Edges Weight: %.2f, True Positive Weight: %.2f",
          weightActualPositive, weightPredictedPositive, weightTruePositive
      );
      line4 = String.format(
          "\nWeighted Precision: %.3f, Weighted Recall: %.3f, Weighted F1: %.3f",
          weightedPrecision(), weightedRecall(), weightedF1()
      );
    }
    return line1 + "\n" + line2 + line3 + line4;
  }
}
