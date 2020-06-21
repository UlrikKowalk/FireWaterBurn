//import java.io.BufferedWriter;
//import java.io.File;
//import java.io.FileWriter;
//import java.util.ArrayList;

class Localisation {

    private int sensors;
    private int blocksize;
    //private int hopsize;
    //private double[][] coordinates;
    //private static final double SPEED_OF_SOUND = 343;
    //private static final int MAXSOURCES = 5;
    private Complex[][][] psd;
    //private int nBins;
    //private Complex alphaCpx;
    //private Complex oneMinusAlphaCpx;
    //private int samplerate;
    private Complex[][][] mTheta;
    //private double[] frequencies;
    //private Complex[][] complexZero;
    //private Complex[] complexSensorVector;
    //private Complex complexZeroNum;
    private double alpha;
    private double oneMinusAlpha;
    //private double overlap;
    private int nLowerBin, nUpperBin;
    private int blocksize_half;
    private int num_theta;
    private double[] v_P_abs_sum;
    //private double[] angles;
    private boolean firstBlock = true;
    //private double nMinMag;
    //private double dt;
    private double[] theta;
    //private int max_theta;
    //private int step_theta;
    private SteeringVector steeringVector;
    

    //private int bTMPBLOCK = 0;

    //private SourceManager sourceManager;

    

    public Localisation(int sensors, int samplerate, int blocksize) {

        this.sensors = sensors;
        //this.samplerate = samplerate;
        this.blocksize = blocksize;

        this.nLowerBin = 100;
        this.blocksize_half = (int) (this.blocksize / 2 + 1);
        this.nUpperBin = this.blocksize_half;
        
        this.alpha = 0.4f;
        //this.nMinMag = 1.2f;
        this.oneMinusAlpha = 1.0f - this.alpha;
        //this.alphaCpx = new Complex(this.alpha, 0);
        //this.oneMinusAlphaCpx = new Complex((1 - this.alpha), 0);

        double max_theta = 360;
        double step_theta = 10;
        this.num_theta = (int)(max_theta/step_theta);

        this.theta = new double[this.num_theta];
        int idxTheta = 0;
        int idx = 0;
        while (idxTheta < max_theta) {
            this.theta[idx] = 2 * Math.PI * idxTheta / 360f;
            idx++;
            idxTheta += step_theta;
        }

        //this.hopsize = (int) (this.blocksize * (1.0 - this.overlap));
        this.steeringVector = new SteeringVector(sensors, samplerate, blocksize);
        this.mTheta = this.steeringVector.generateDelayTensor(this.theta);

        this.v_P_abs_sum = new double[this.num_theta];

        

        System.out.println("Theta: " + this.theta.length);

    }

    public int getNumTheta() {
        return this.num_theta;
    }


    /*public double[][] runExperiment(double[][] audioData) {

        

        

        


     


        
        //this.coordinates = establishCoordinates();

        

        

        //ResultWriter resultWriter = new ResultWriter("threshold.txt");
        //for (int iBlock = 0; iBlock < result.length; iBlock++) {
        //    resultWriter.write(result[iBlock]);
        //}

        return result;
    }*/


    private Complex[][] getSteeringVector(int iTheta) {
        return this.mTheta[iTheta];
    }



    

    private Complex getRealTrace(Complex[][][] matrix, int iFreq) {

        double real = 0f;

        for (int iRow = 0; iRow < this.sensors; iRow++) {
            real += matrix[iRow][iRow][iFreq].re;
        }
        return new Complex(real, 0f);
    }

    /*private double[][] generateEye(int size) {
        double[][] eye = new double[size][size];
        for (int iSize = 0; iSize < size; iSize++) {
            eye[iSize][iSize] = 1;
        }
        return eye;
    }*/

    /*private void generateComplexZero(int size) {
        this.complexZero = new Complex[size][size];
        Complex zero = new Complex(0,0);
        for (int iRow = 0; iRow < size; iRow++) {
            for (int iCol = 0; iCol < size; iCol++) {
                this.complexZero[iRow][iCol] = zero;
            }
        }
    }*/

    private Complex[][] generateEyeTimesComplex(int size, Complex number) {
        Complex[][] eye = new Complex[size][size];
        for (int iRow = 0; iRow < size; iRow++) {
            for (int iCol = 0; iCol < size; iCol++) {
                eye[iRow][iCol] = new Complex(0,0);
            }
        }
        for (int iRow = 0; iRow < size; iRow++) {
            eye[iRow][iRow] = number;
        }
        return eye;
    }

    /*private int argMax(double[] inVec) {
        int arg = 0;
        double max = inVec[0];
        for (int idx = 0; idx < inVec.length; idx++) {
            if (inVec[idx] > max) {
                max = inVec[idx];
                arg = idx;
            }
        }
        return arg;
    }*/

    /*private double max(double[] inVec) {
        double max = inVec[0];
        for (int idx = 0; idx < inVec.length; idx++) {
            if (inVec[idx] > max) {
                max = inVec[idx];
            }
        }
        return max;
    }*/

    /*private void generateComplexSensorVector() {
        this.complexSensorVector = new Complex[this.sensors];
        for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
            this.complexSensorVector[iSensor] = new Complex(0,0);
        }
    }*/

    /*protected void writeResult(String data, String filename) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename, true));
            writer.append(data);
            writer.close();
        } catch (Exception e) { }
    }*/




    /*private double[][] estimateDoa(double[][] sensor_data) {

        

        
        
        
        

        //this.complexZeroNum = new Complex(0,0);

        //generateComplexZero(this.sensors);
        
        //generateComplexSensorVector();

        

        
        




            double[] directions = sourceDirections(spec);



            //for (int iPeak = 0; iPeak < this.num_theta + 10 + 2*MAXSOURCES; iPeak++) {
            //    mAbs_sum[iBlock][iPeak] = directions[iPeak];
            //}



            
        }

        
    return mAbs_sum;
}*/


    public double[][] sourceDirections(Complex[][] spec) {

        // Generation of the PSD

        if (this.firstBlock) {

            this.psd = new Complex[this.sensors][this.sensors][blocksize_half];

            for (int iFreq = this.nLowerBin; iFreq < this.nUpperBin; iFreq++) {
                for (int iRow = 0; iRow < this.sensors; iRow++) {
                    for (int iCol = 0; iCol < this.sensors; iCol++) {
                        this.psd[iRow][iCol][iFreq] = new Complex(
                            spec[iRow][iFreq].re*spec[iCol][iFreq].re + spec[iRow][iFreq].im*spec[iCol][iFreq].im,
                            spec[iRow][iFreq].im*spec[iCol][iFreq].re - spec[iRow][iFreq].re*spec[iCol][iFreq].im);
                    }
                }
            }

        } else {

            for (int iFreq = this.nLowerBin; iFreq < this.nUpperBin; iFreq++) {
                for (int iRow = 0; iRow < this.sensors; iRow++) {
                    for (int iCol = 0; iCol < this.sensors; iCol++) {

                        this.psd[iRow][iCol][iFreq].re = this.oneMinusAlpha * 
                            (spec[iRow][iFreq].re*spec[iCol][iFreq].re + spec[iRow][iFreq].im*spec[iCol][iFreq].im) +
                            this.alpha * this.psd[iRow][iCol][iFreq].re;
                        this.psd[iRow][iCol][iFreq].im = this.oneMinusAlpha * 
                            (spec[iRow][iFreq].im*spec[iCol][iFreq].re - spec[iRow][iFreq].re*spec[iCol][iFreq].im) +
                            this.alpha * this.psd[iRow][iCol][iFreq].im;        
                    }
                }
            }
        }


        // Testing of various azimuth angles

        double nP = 0;

        for (int iTheta = 0; iTheta < num_theta; iTheta++) {

            nP = 0;

            Complex[][] vectorA = getSteeringVector(iTheta);

            // Power formula 27
            for (int iBin = this.nLowerBin; iBin < this.nUpperBin; iBin++) {

                Complex trace = getRealTrace(this.psd, iBin);
            
                Complex[][] eyeTimesTraceMinusPSD = generateEyeTimesComplex(this.sensors, trace);

                for (int iRow = 0; iRow < this.sensors; iRow++) {
                    for (int iCol = 0; iCol < this.sensors; iCol++) {

                        eyeTimesTraceMinusPSD[iRow][iCol].re -= this.psd[iRow][iCol][iBin].re;
                        eyeTimesTraceMinusPSD[iRow][iCol].im -= this.psd[iRow][iCol][iBin].im;
                    }
                }

                Complex[] tmp = new Complex[this.sensors];
                for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                    tmp[iSensor] = new Complex(0,0);
                }

                for (int iRow = 0; iRow < this.sensors; iRow++) {
                    for (int iCol = 0; iCol < this.sensors; iCol++) {

                        tmp[iRow].re += (eyeTimesTraceMinusPSD[iRow][iCol].re * vectorA[iCol][iBin].re - 
                                            eyeTimesTraceMinusPSD[iRow][iCol].im * vectorA[iCol][iBin].im); 
                        //tmp[iRow].im += (eyeTimesTraceMinusPSD[iRow][iCol].im * vectorA[iCol][iBin].re + 
                        //                    eyeTimesTraceMinusPSD[iRow][iCol].re * vectorA[iCol][iBin].im); 
                    }
                }

                Complex tmp2 = new Complex(0, 0);
                for (int iRow = 0; iRow < this.sensors; iRow++) {
                    for (int iCol = 0; iCol < this.sensors; iCol++) {

                        tmp2.re += vectorA[iCol][iBin].re*tmp[iRow].re + 
                            vectorA[iCol][iBin].im*tmp[iRow].im;
                        //tmp2.im += vectorA[iCol][iBin].re*tmp[iRow].im - 
                        //    vectorA[iCol][iBin].im*tmp[iRow].re; 

                    }
                }

                // Power summation over all (interesting) bins
                //nP += tmp2.abs();
                nP += tmp2.re;
            }
            
            v_P_abs_sum[iTheta] = -1f / nP;

        }

        // Loudness Normalization
        double rms = Loudness.rms(v_P_abs_sum);
        for (int iTheta = 0; iTheta < this.num_theta; iTheta++) {
            v_P_abs_sum[iTheta] /= rms;
        }

        // Peak Picking
        int[] vPeaks = FindPeaks.findPeaks(v_P_abs_sum);

        // Quadratic Interpolation of Peak Positions
        double[][] mRealPeaks = QuadraticInterpolation.findRealPeaks(v_P_abs_sum, vPeaks);


        
        


        // Write all data to text file
        /*double[] tmp = new double[this.num_theta + 10 + 2 * SourceManager.MAX_SOURCES];

        // Regular Source Distribution
        for (int iTheta = 0; iTheta < this.num_theta; iTheta++) {
            tmp[iTheta] = v_P_abs_sum[iTheta];
        }
        // Interpolated Peak Results
        for (int iTheta = 0; iTheta < mRealPeaks[0].length; iTheta++) {
            if (mRealPeaks[iTheta][1] > this.nMinMag) {
                tmp[this.num_theta + iTheta] = mRealPeaks[iTheta][0];
            }
        }
        // Kalman Filtered Sources
        for (int iSource = 0; iSource < 2*SourceManager.MAX_SOURCES; iSource++) {
            tmp[this.num_theta + 10 + iSource] = vKalmanPeaks[iSource];
        }*/
        if (this.firstBlock) {
            this.firstBlock = false;
        }
        
        return mRealPeaks;
    }



    /*private double calculateNorm(double[] vector) {
        return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    }*/

    /*private double NaNtoZero(double number) {
        if (Double.isNaN(number)) {
            return 0;
        } else {
            return number;
        }
    }*/


    /*private double vectorAngle(double[] unit, double[] sensor) {

        // Calculate the Angle between two Vectors

        double dotProduct = unit[0] * sensor[0] + unit[1] * sensor[1] + unit[2] * sensor[2];
        double normUnit = calculateNorm(unit);
        double normSensor = calculateNorm(sensor);
        double angle = Math.acos((double) (dotProduct / (normUnit * normSensor)));

        if (sensor[1] < unit[1]) {
            angle = -angle;
        }

        return angle;
    }*/

    /*private double[][] establishCoordinates() {

        // Generate Coordinates for MatrixVoice Array

        double[][] coordinates = new double[this.sensors][3];
        double omega;

        // Not necessary, just feels right
        coordinates[0][0] = 0;
        coordinates[0][1] = 0;
        coordinates[0][2] = 0;

        for (int iSensor = 1; iSensor < this.sensors; iSensor++) {
            omega = -2f * Math.PI * (iSensor - 1) / (this.sensors - 1) - Math.PI;
            coordinates[iSensor][0] = Math.PI * Math.cos(omega) * 0.01;
            coordinates[iSensor][1] = Math.PI * Math.sin(omega) * 0.01;
            coordinates[iSensor][2] = 0;
        }

        return coordinates;
    }*/
    
    /*private double[] calculateDelays(double theta) {

        //Calculate the time delay at each microphone for a given angle theta
        
        double[][] coords = new double[this.sensors][3];

        // Center coordinates around first entry
        for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
            coords[iSensor][0] = this.coordinates[iSensor][0] - this.coordinates[this.sensors - 1][0];
            coords[iSensor][1] = this.coordinates[iSensor][1] - this.coordinates[this.sensors - 1][1];
            coords[iSensor][2] = this.coordinates[iSensor][2] - this.coordinates[this.sensors - 1][2];
        }

        double[] v_unit = new double[] { 1, 0, 0 };
        double[] delays = new double[sensors];

        for (int iSensor = 0; iSensor < sensors; iSensor++) {

            delays[iSensor] = NaNtoZero(calculateNorm(coords[iSensor])
                    * Math.cos(theta + vectorAngle(v_unit, coords[iSensor])) / this.SPEED_OF_SOUND);
        }

        return delays;
    }*/

    /*private void generateDelayTensor(double[] theta, int blocklen_half) {
        
        this.mTheta = new Complex[theta.length][this.sensors][blocklen_half];
        double[] v_tau_theta = new double[this.sensors];

        for (int iTheta = 0; iTheta < theta.length; iTheta++) {
            // Calculation of Delays
            v_tau_theta = calculateDelays(theta[iTheta]);
            
            // Creation of Steering Vector
            for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                for (int iFreq = 0; iFreq < blocklen_half; iFreq++) {
                    this.mTheta[iTheta][iSensor][iFreq] = new Complex(0,0);
                    this.mTheta[iTheta][iSensor][iFreq].re = Math.cos(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                    this.mTheta[iTheta][iSensor][iFreq].im = Math.sin(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                }
            }
        }
    }*/


}