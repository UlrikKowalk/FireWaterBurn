package com.jadehs.viwer;

public class Source {

    private KalmanNew kalman;
    private double nAngle;
    private double nSpread;
    private double acceptance_probability = 8.0f; //4.0f
    private int blocksWithoutUpdate = 0;
    private double spatial_smoothing_constant = 0.85f; //0.85f
    private boolean isNew;

    public Source(double angle, double dt, double nSpread) {

        this.nAngle = angle;
        this.nSpread = nSpread;
        this.kalman = new KalmanNew(dt, this.nAngle);
        this.isNew = true;

    }


    public double iterate(double nCandidate) {
        this.isNew = false;
        this.blocksWithoutUpdate = 0;
        this.nAngle = (1.0f - this.spatial_smoothing_constant) * this.kalman.iterate(nCandidate) + this.spatial_smoothing_constant * this.nAngle;
        return this.nAngle;
    }

    public double iterateWeighted(double[] vCandidates) {

        // Account for cyclic range 
        int nCandidates = vCandidates.length;
                
        double[] vTempCandidates = new double[3*nCandidates];
        for (int iCandidate = 0; iCandidate < nCandidates; iCandidate++) {
            vTempCandidates[iCandidate] = vCandidates[iCandidate] - this.nSpread;
            vTempCandidates[iCandidate + nCandidates] = vCandidates[iCandidate];
            vTempCandidates[iCandidate + 2*nCandidates] = vCandidates[iCandidate] + this.nSpread;
        } 

        double[] vWeights = new double[3*nCandidates];
        double nSumWeights = 0.0f;
        double tmp = 0.0f;
        for (int iWeight = 0; iWeight < 3*nCandidates; iWeight++) {
            tmp = pdf(vTempCandidates[iWeight], this.nAngle, this.acceptance_probability);
            vWeights[iWeight] = tmp;
            nSumWeights += tmp;
        }

        double nEstiPos = 0.0f;
        for (int iCandidate = 0; iCandidate < nCandidates; iCandidate++) {
            nEstiPos += vCandidates[iCandidate] * vWeights[iCandidate];
        }
        for (int iCandidate = 0; iCandidate < nCandidates; iCandidate++) {
            nEstiPos += vCandidates[iCandidate] * vWeights[iCandidate + nCandidates];
        }
        for (int iCandidate = 0; iCandidate < nCandidates; iCandidate++) {
            nEstiPos += vCandidates[iCandidate] * vWeights[iCandidate + 2*nCandidates];
        }
        nEstiPos /= nSumWeights;
        
        return iterate(nEstiPos);
    }

 
    public double noUpdate() {
        this.blocksWithoutUpdate++;
        return this.kalman.getEstimate();
    }

    public boolean getIsNew() {
        return this.isNew;
    }

    public int getBlocksWithoutUpdate() {
        return this.blocksWithoutUpdate;
    }

    public double getAngle() {
        return this.kalman.getEstimate();
    }

    public double getPdf(double x) {
        return pdf(x, getAngle(), this.acceptance_probability);
    }

    // return pdf(x) = standard Gaussian pdf
    public static double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return pdf(x, mu, sigma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
        return pdf((x - mu) / sigma) / sigma;
    }



}