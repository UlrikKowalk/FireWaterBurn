package com.JadeHS.Viwer;

class SteeringVector {

    private static final double SPEED_OF_SOUND = 343;
    private int sensors;
    private int blocksize;
    private int blocklen_half;
    private double[] frequencies;
    private double[][] coordinates;

    public SteeringVector(int sensors, int samplerate, int blocksize) {

        this.sensors = sensors;
        this.blocksize = blocksize;
        this.blocklen_half = (int) (this.blocksize / 2 + 1);

        double samplerateDividedByBlocksize = (double) samplerate / blocksize;
        
        this.frequencies = new double[blocklen_half];
        for (int iFreq = 0; iFreq < blocklen_half; iFreq++) {
            this.frequencies[iFreq] = iFreq * samplerateDividedByBlocksize;
        }
        
        establishCoordinates();
    }

    public double[][] getCoordinates() {
        return this.coordinates;
    }

    private void establishCoordinates() {

        //Generate Coordinates for MatrixVoice Array
        
        this.coordinates = new double[this.sensors][3];
        double omega;

        // Not necessary, just feels right
        this.coordinates[0][0] = 0;
        this.coordinates[0][1] = 0;
        this.coordinates[0][2] = 0;

        for (int iSensor = 1; iSensor < this.sensors; iSensor++) {
            omega = -2f * Math.PI * (iSensor - 1) / (this.sensors - 1) - Math.PI;
            this.coordinates[iSensor][0] = Math.PI * Math.cos(omega) * 0.01f;
            this.coordinates[iSensor][1] = Math.PI * Math.sin(omega) * 0.01f;
            this.coordinates[iSensor][2] = 0f;
        }
    }

    public double[] calculateDelays(double theta) {

        // Calculate the time delay at each microphone for a given angle theta

        double[][] coords = new double[this.sensors][3];

        // Center coordinates around first entry
        for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
            coords[iSensor][0] = this.coordinates[iSensor][0] - this.coordinates[this.sensors - 1][0];
            coords[iSensor][1] = this.coordinates[iSensor][1] - this.coordinates[this.sensors - 1][1];
            coords[iSensor][2] = this.coordinates[iSensor][2] - this.coordinates[this.sensors - 1][2];
        }

        double[] v_unit = new double[] { 1f, 0f, 0f };
        double[] delays = new double[sensors];

        for (int iSensor = 0; iSensor < sensors; iSensor++) {

            delays[iSensor] = NaNtoZero(calculateNorm(coords[iSensor])
                    * Math.cos(theta + vectorAngle(v_unit, coords[iSensor])) / SteeringVector.SPEED_OF_SOUND); 
                           
        }

        return delays;
    }

    public double[][][] generateDelayTensor(double theta) {
        
        // Sensors, Bins, Real/Imag
        double[][][] mTheta = new double[this.sensors][blocklen_half][2];
        double[] v_tau_theta = new double[this.sensors];

        // Calculation of Delays
        v_tau_theta = calculateDelays(theta);
        
        // Creation of Steering Vector
        for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
            for (int iFreq = 0; iFreq < blocklen_half; iFreq++) {
                mTheta[iSensor][iFreq][0] = Math.cos(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                mTheta[iSensor][iFreq][1] = Math.sin(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
            }
        }
        
        return mTheta;
    }

    public Complex[][][] generateDelayTensor(double[] theta) {
        
        Complex[][][] mTheta = new Complex[theta.length][this.sensors][blocklen_half];
        double[] v_tau_theta = new double[this.sensors];

        for (int iTheta = 0; iTheta < theta.length; iTheta++) {
            // Calculation of Delays
            v_tau_theta = calculateDelays(theta[iTheta]);
            
            // Creation of Steering Vector
            for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                for (int iFreq = 0; iFreq < blocklen_half; iFreq++) {
                    mTheta[iTheta][iSensor][iFreq] = new Complex(0,0);
                    mTheta[iTheta][iSensor][iFreq].re = Math.cos(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                    mTheta[iTheta][iSensor][iFreq].im = Math.sin(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                }
            }
        }
        return mTheta;
    }

    private double vectorAngle(double[] unit, double[] sensor) {

        // Calculate the Angle between two Vectors

        double dotProduct = unit[0] * sensor[0] + unit[1] * sensor[1] + unit[2] * sensor[2];
        double normUnit = calculateNorm(unit);
        double normSensor = calculateNorm(sensor);
        double angle = Math.acos((double) (dotProduct / (normUnit * normSensor)));

        if (sensor[1] < unit[1]) {
            angle = -angle;
        }

        return angle;
    }

    private double calculateNorm(double[] vector) {
        return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    }

    private double NaNtoZero(double number) {
        if (Double.isNaN(number)) {
            return 0;
        } else {
            return number;
        }
    }
}