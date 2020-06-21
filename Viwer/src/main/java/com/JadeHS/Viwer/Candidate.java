package com.jadehs.viwer;

public class Candidate {

    private double theta = -255f;
    private double probability = -255f;

    public Candidate(double theta, double probability) {
        this.theta = theta;
        this.probability = probability;
    }

    public double getTheta() {
        return this.theta;
    }

    public double getProbability() {
        return this.probability;
    }

    public double[] getData() {
        double[] data = {this.theta, this.probability};
        return data;
    }
    



}