
import java.util.*; 

class Kalman {

    private static int numSources = 0;
    private final static int NUMBLOCKSTILLDEATH = 10;

    private int numBlocksRemaining;

    private boolean bFirstBlock = true;
    private double dt;
    private double nObservationError = 10*1e-3f;
    private double nObservation;

    private double[] vKalmanGain = {0.0f, 0.0f};
    private double[] vEstimate = {0.0f, 0.5f};
    private double[] vH = {1.0f, 0.0f};

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
        this.numBlocksRemaining = NUMBLOCKSTILLDEATH;
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

        /*if (this.bFirstBlock) {

            this.bFirstBlock = false;
            this.vEstimate[0] = nObservation;
            return nObservation;

        } else {*/

            // State Covariance Update
            updateStateCovariance();

            // Kalman Gain
            calculateKalmanGain();

            // Process Error Update
            updateProcessError();

            // Estimation
            return estimate();
        
        //}
    }

    private void makeObservation(double nNewMeasurement) {
        this.nObservation = nNewMeasurement + this.nObservationError * oRandom.nextGaussian();
    }

    private void updateStateCovariance() {

        /*this.mQ[0][0] += oRandom.nextGaussian();
        this.mQ[0][1] += oRandom.nextGaussian();
        this.mQ[1][0] += oRandom.nextGaussian();
        this.mQ[1][1] += oRandom.nextGaussian();*/

        this.mStateCovariance = addTwoByTwo(multiplyTwoByTwoTransposed(multiplyTwoByTwo(this.mA, this.mProcessError), this.mA), this.mQ);
    }

    private void calculateKalmanGain() {
        double tmp = (this.mStateCovariance[0][0] + this.nObservationError * this.nObservationError);
        this.vKalmanGain[0] = this.mStateCovariance[0][0] / tmp;
        this.vKalmanGain[1] =  this.mStateCovariance[1][0] / tmp;

        double a = 1.0f;
    }

    private void updateProcessError() {
        this.mProcessError[0][0] = this.mStateCovariance[0][0] - this.vKalmanGain[0] * this.mStateCovariance[0][0];
        this.mProcessError[0][1] = this.mStateCovariance[0][1] - this.vKalmanGain[0] * this.mStateCovariance[0][1];
        this.mProcessError[1][0] = this.mStateCovariance[1][0] - this.vKalmanGain[1] * this.mStateCovariance[0][0];
        this.mProcessError[1][1] = this.mStateCovariance[1][1] - this.vKalmanGain[1] * this.mStateCovariance[0][1];
    }

    private double estimate() {

        double[] tmp = {this.mA[0][0] * this.vEstimate[0] + this.mA[0][1] * this.vEstimate[1], this.mA[1][0] * this.vEstimate[0] + this.mA[1][1] * this.vEstimate[1]};

        //double tmpDifference = this.nObservation - this.vEstimate[0] * (this.mA[0][0] * this.vEstimate[0] + this.mA[0][1] * this.vEstimate[1]); 
        double tmpDifference = this.nObservation - (this.mA[0][0] * this.vEstimate[0] + this.mA[0][1] * this.vEstimate[1]); 


        this.vEstimate[0] = tmp[0] + this.vKalmanGain[0] * tmpDifference;
        this.vEstimate[1] = tmp[1] + this.vKalmanGain[1] * tmpDifference;

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


}