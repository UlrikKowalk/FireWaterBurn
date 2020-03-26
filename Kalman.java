
import java.util.Random;

class Kalman {

    private static int numSources = 0;
    private final static int NUMSECONDSTILLDEATH = 5;

    private int numBlocksRemaining;

    private double dt;
    private double nObservationError = 1e-3f;
    private double nObservation;

    private double[] vKalmanGain = {0.0f, 0.0f};
    private double[] vEstimate = {0.0f, 0.0f};

    private double[][] mA;
    private double[][] mQ = newZeroMatrix(2, 2);
    private double[][] mProcessError = {{1.0f, 0.0f}, {0.0f, 1.0f}};
    private double[][] mStateCovariance = newZeroMatrix(2, 2);

    private Random oRandom = new Random();


    public Kalman(double dt) {

        numSources++;
        this.dt = dt;
        this.mA = new double[][] {{1.0f, this.dt}, {0.0f, 1.0f}};
        this.resetRemaining();

    }

    public static int getNumberOfSources() {
        return numSources;
    }

    public static double calculateDt(int blocksize, int samplingrate) {
        return ((double) blocksize) / samplingrate; 
    }

    private void resetRemaining() {
        this.numBlocksRemaining = (int) (NUMSECONDSTILLDEATH / this.dt);
    }

    private void decreaseRemaining() {
        this.numBlocksRemaining--;
    }

    public int getNumBlocksRemaining() {
        return this.numBlocksRemaining;
    }

    public double iterate(double nNewMeasurement) {

        // Observed State
        makeObservation(nNewMeasurement);

        // State Covariance Update
        updateStateCovariance();

        // Kalman Gain
        calculateKalmanGain();

        // Process Error Update
        updateProcessError();

        // Estimation
        return estimate();

    }

    public double iterateWeighted(double[] vCandidates, double[] vTheta, double nWidth) {

        // Return a weighted mean source candidate

        // Account for cyclic range 
        int nCandidates = vCandidates.length;
        double nSpread = Math.abs(vTheta[vTheta.length-1] - vTheta[0]);
        double[] vTempCandidates = new double[3*nCandidates];
        for (int iCandidate = 0; iCandidate < nCandidates; iCandidate++) {
            vTempCandidates[iCandidate] = vCandidates[iCandidate] - nSpread;
            vTempCandidates[iCandidate + nCandidates] = vCandidates[iCandidate];
            vTempCandidates[iCandidate + 2*nCandidates] = vCandidates[iCandidate] + nSpread;
        } 

        double[] vWeights = new double[3*nCandidates];
        double nSumWeights = 0.0f;
        double tmp = 0.0f;
        for (int iWeight = 0; iWeight < 3*nCandidates; iWeight++) {
            tmp = pdf(vTempCandidates[iWeight], this.vEstimate[0], nWidth * (5 * this.vKalmanGain[0] + 0.5f));
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

        return this.iterate(nEstiPos);

    }

    private void makeObservation(double nNewMeasurement) {

        this.nObservation = nNewMeasurement + this.nObservationError * oRandom.nextGaussian();
        
    }

    private void updateStateCovariance() {

        this.mStateCovariance = addTwoByTwo(multiplyTwoByTwoTransposed(multiplyTwoByTwo(this.mA, this.mProcessError), this.mA), this.mQ);

    }

    private void calculateKalmanGain() {

        double tmp = (this.mStateCovariance[0][0] + this.nObservationError * this.nObservationError);
        this.vKalmanGain[0] = this.mStateCovariance[0][0] / tmp;
        this.vKalmanGain[1] =  this.mStateCovariance[1][0] / tmp;

    }

    private void updateProcessError() {

        this.mProcessError[0][0] = this.mStateCovariance[0][0] - this.vKalmanGain[0] * this.mStateCovariance[0][0];
        this.mProcessError[0][1] = this.mStateCovariance[0][1] - this.vKalmanGain[0] * this.mStateCovariance[0][1];
        this.mProcessError[1][0] = this.mStateCovariance[1][0] - this.vKalmanGain[1] * this.mStateCovariance[0][0];
        this.mProcessError[1][1] = this.mStateCovariance[1][1] - this.vKalmanGain[1] * this.mStateCovariance[0][1];

    }

    private double estimate() {

        double[] tmp = {this.mA[0][0] * this.vEstimate[0] + this.mA[0][1] * this.vEstimate[1], this.mA[1][0] * this.vEstimate[0] + this.mA[1][1] * this.vEstimate[1]};

        double tmpDifference = this.nObservation - (this.mA[0][0] * this.vEstimate[0] + this.mA[0][1] * this.vEstimate[1]); 

        this.vEstimate[0] = tmp[0] + this.vKalmanGain[0] * tmpDifference;
        this.vEstimate[1] = tmp[1] + this.vKalmanGain[1] * tmpDifference;

        return this.vEstimate[0];

    }

    public double getEstimate() {
        return this.vEstimate[0];
    }

    private double[][] newZeroMatrix(int nRows, int nCols) {
        double[][] mZero = new double[nRows][nCols];
        for (int iCol = 0; iCol < nCols; iCol++) {
            for (int iRow = 0; iRow < nRows; iRow++) {
                mZero[iRow][iCol] = 0.0f;
            }
        }
        return mZero;
    }

    private double[][] addTwoByTwo(double[][] mA, double[][] mB) {

        // 2x2 Matrix Addition 

        double[][] mResult = {{0.0f, 0.0f}, {0.0f, 0.0f}};
        
        mResult[0][0] = mA[0][0] + mB[0][0];
        mResult[0][1] = mA[0][1] + mB[0][1];
        mResult[1][0] = mA[1][0] + mB[1][0];
        mResult[1][1] = mA[1][1] + mB[1][1];

        return mResult;
    }

    private double[][] multiplyTwoByTwo(double[][] mA, double[][] mB) {

        // 2x2 Matrix Multiplication 

        double[][] mResult = {{0.0f, 0.0f}, {0.0f, 0.0f}};
        
        mResult[0][0] = mA[0][0]*mB[0][0] + mA[0][1]*mB[1][0];
        mResult[0][1] = mA[0][0]*mB[0][1] + mA[0][1]*mB[1][1];
        mResult[1][0] = mA[1][0]*mB[0][0] + mA[1][1]*mB[1][0];
        mResult[1][1] = mA[1][0]*mB[0][1] + mA[1][1]*mB[1][1];

        return mResult;
    }

    private double[][] multiplyTwoByTwoTransposed(double[][] mA, double[][] mB) {

        // 2x2 Matrix Multiplication 

        double[][] mResult = {{0.0f, 0.0f}, {0.0f, 0.0f}};
        
        mResult[0][0] = mA[0][0]*mB[0][0] + mA[0][1]*mB[0][1];
        mResult[0][1] = mA[0][0]*mB[1][0] + mA[0][1]*mB[1][1];
        mResult[1][0] = mA[1][0]*mB[0][0] + mA[1][1]*mB[0][1];
        mResult[1][1] = mA[1][0]*mB[1][0] + mA[1][1]*mB[1][1];

        return mResult;
    }

     // return pdf(x) = standard Gaussian pdf
     public static double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
        return pdf((x - mu) / sigma) / sigma;
    }

}