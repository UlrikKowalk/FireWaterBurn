import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

class Localisation {

    private int sensors;
    private double[][] coordinates;
    private double SPEED_OF_SOUND = 343;
    private Complex[][][] psd;
    private Complex[][] spec;
    private int nBins;
    private Complex alphaCpx;
    private Complex oneMinusAlphaCpx;
    private int samplerate;
    private Complex[][][] mTheta;
    private double[] frequencies;
    private Complex[][] complexZero;
    private Complex[] complexSensorVector;
    private Complex complexZeroNum;
    private double alpha;
    private double oneMinusAlpha;
    private int nLowerBin, nUpperBin;
    private int blocklen_half;

    public Localisation() {

    }

    public double[] performExperiment(double[][] audioData, int samplerate, int blocksize, double overlap) {

        this.sensors = audioData.length;
        this.samplerate = samplerate;

        this.nLowerBin = 0;
        this.nUpperBin = (int) (blocksize / 2 + 1);
        this.blocklen_half = (int) (blocksize / 2 + 1);
        this.alpha = 0.9f;
        this.oneMinusAlpha = 1.0f - this.alpha;
        this.alphaCpx = new Complex(this.alpha, 0);
        this.oneMinusAlphaCpx = new Complex((1 - this.alpha), 0);

        this.coordinates = establishCoordinates();

        double[] result = estimateDoa(audioData, samplerate, coordinates, blocksize, overlap);

        ResultWriter resultWriter = new ResultWriter("threshold.txt");
        for (int iBlock = 0; iBlock < result.length; iBlock++) {
            resultWriter.write(result[iBlock]);
        }

        return result;

        // # todo: check shape
        // for azimuth in all_azimuth:
        // sensor_data = self.generate_room(azimuth, coordinates, distance, raw_data,
        // samplerate)
        // estimated_doa.append(np.median(self.estimate_doa(sensor_data, samplerate,
        // coordinates, blocksize, overlap)))

    }

    private double vectorAngle(double[] unit, double[] sensor) {

        /**
         * Calculate the Angle between two Vectors
         */

        double dotProduct = unit[0] * sensor[0] + unit[1] * sensor[1] + unit[2] * sensor[2];
        double normUnit = calculateNorm(unit);
        double normSensor = calculateNorm(sensor);
        double angle = Math.acos((double) (dotProduct / (normUnit * normSensor)));

        if (sensor[1] < unit[1]) {
            angle = -angle;
        }

        return angle;
    }

    private double[][] establishCoordinates() {

        /**
         * Generate Coordinates for MatrixVoice Array
         */

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
    }
    
    private double[] calculateDelays(double theta) {

        /**
         * Calculate the time delay at each microphone for a given angle theta
         */

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
    }

    private void generateDelayTensor(double[] theta, int blocklen_half) {
        
        this.mTheta = new Complex[theta.length][this.sensors][blocklen_half];
        double[] v_tau_theta = new double[this.sensors];

        for (int iTheta = 0; iTheta < theta.length; iTheta++) {
            // Calculation of Delays
            v_tau_theta = calculateDelays(theta[iTheta]);
            
            // Creation of Steering Vector
            for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                for (int iFreq = 0; iFreq < blocklen_half; iFreq++) {
                    //mTheta[iTheta][iSensor][iFreq] = (new Complex(0,
                    //        -2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor])).exp();

                    //mTheta[iTheta][iSensor][iFreq] = new Complex(Math.cos(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]), 
                    //    Math.sin(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]));
                    this.mTheta[iTheta][iSensor][iFreq] = new Complex(0,0);
                    this.mTheta[iTheta][iSensor][iFreq].re = Math.cos(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                    this.mTheta[iTheta][iSensor][iFreq].im = Math.sin(-2f * Math.PI * this.frequencies[iFreq] * v_tau_theta[iSensor]);
                }
            }
        }
    }

    private Complex[][] getSteeringVector(int iTheta) {
        return this.mTheta[iTheta];
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

    private double[] hann(int blocksize) {
        double[] window = new double[blocksize];
        for (int iSample = 0; iSample < blocksize; iSample++) {
            window[iSample] = 0.5 * (1 - Math.cos(2 * Math.PI * iSample / (blocksize - 1)));
        }
        return window;
    }

    private Complex getRealTrace(Complex[][][] matrix, int iFreq) {

        double real = 0f;
        //double imag = 0f;

        for (int iRow = 0; iRow < this.sensors; iRow++) {
            //trace = trace.plus(matrix[iRow][iRow][iFreq]);
            real += matrix[iRow][iRow][iFreq].re;
            //imag += matrix[iRow][iRow][iFreq].im;
        }
        //return new Complex(real, imag);
        return new Complex(real, 0f);
    }

    private double[][] generateEye(int size) {
        double[][] eye = new double[size][size];
        for (int iSize = 0; iSize < size; iSize++) {
            eye[iSize][iSize] = 1;
        }
        return eye;
    }

    private void generateComplexZero(int size) {
        this.complexZero = new Complex[size][size];
        Complex zero = new Complex(0,0);
        for (int iRow = 0; iRow < size; iRow++) {
            for (int iCol = 0; iCol < size; iCol++) {
                this.complexZero[iRow][iCol] = zero;
            }
        }
    }

    private Complex[][] generateEyeTimesComplex(int size, Complex number) {
        //Complex[][] eye = this.complexZero.clone();
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

    private int argMax(double[] inVec) {
        int arg = 0;
        double max = inVec[0];
        for (int idx = 0; idx < inVec.length; idx++) {
            if (inVec[idx] > max) {
                max = inVec[idx];
                arg = idx;
            }
        }
        return arg;
    }

    private double max(double[] inVec) {
        double max = inVec[0];
        for (int idx = 0; idx < inVec.length; idx++) {
            if (inVec[idx] > max) {
                max = inVec[idx];
            }
        }
        return max;
    }

    private void generateComplexSensorVector() {
        this.complexSensorVector = new Complex[this.sensors];
        for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
            complexSensorVector[iSensor] = new Complex(0,0);
        }
    }

    protected void writeResult(String data, String filename) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename, true));
            writer.append(data);
            writer.close();
        } catch (Exception e) { }
    }

    private double[] estimateDoa(double[][] sensor_data, int samplerate, double[][] coordinates, int blocksize,
            double overlap) {

        int length = sensor_data[0].length;
        double hopsize = blocksize * (1.0 - overlap);
        int num_blocks = (int) Math.floor((length - blocksize) / hopsize);
        double[] window = hann(blocksize);

        int max_theta = 180;
        int step_theta = 10;
        int num_theta = (int)(max_theta/step_theta);
        double theta[] = new double[num_theta];
        int idxTheta = 0;
        int idx = 0;
        while (idxTheta < max_theta) {
            theta[idx] = 2 * Math.PI * idxTheta / 360f;
            idx++;
            idxTheta += step_theta;
        }
        
        double samplerateDividedByBlocksize = (double) samplerate / blocksize;
        this.frequencies = new double[this.blocklen_half];
        for (int iFreq = 0; iFreq < this.blocklen_half; iFreq++) {
            this.frequencies[iFreq] = iFreq * samplerateDividedByBlocksize;
        }

        this.complexZeroNum = new Complex(0,0);

        generateComplexZero(this.sensors);
        generateDelayTensor(theta, this.blocklen_half);
        generateComplexSensorVector();

        this.spec = new Complex[this.sensors][this.blocklen_half];
        this.psd = new Complex[this.sensors][this.sensors][this.blocklen_half];

        double[] angles = new double[num_blocks];
        double[] v_P_abs_sum = new double[num_theta];

        int idx_in, idx_out;
        Complex[] signalBlock = new Complex[blocksize];
        for (int iSample = 0; iSample < blocksize; iSample++) {
            signalBlock[iSample] = new Complex(0,0);
        }

        String filename = "threshold.txt";
        File file = new File(filename);
        file.delete();

        double[][] mAbs_sum = new double[num_blocks][num_theta];

        for (int iBlock = 0; iBlock < num_blocks; iBlock++) {

            idx_in = (int) (iBlock * hopsize);
            idx_out = (int) (idx_in + blocksize);

            // Generation of the Spectrum

            for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                for (int iSample = idx_in; iSample < idx_out; iSample++) {
                    signalBlock[iSample - idx_in].re = sensor_data[iSensor][iSample] * window[iSample - idx_in];
                    signalBlock[iSample - idx_in].im = 0;
                }
                InplaceFFT.fft(signalBlock);
                for (int iFreq = 0; iFreq < this.blocklen_half; iFreq++) {
                    this.spec[iSensor][iFreq] = new Complex(signalBlock[iFreq].re, signalBlock[iFreq].im);
                    //this.spec[iSensor][iFreq] = signalBlock[iFreq];
                }
            }

            // Generation of the PSD

            if (iBlock == 0) {

                this.psd = new Complex[this.sensors][this.sensors][blocklen_half];

                for (int iFreq = this.nLowerBin; iFreq < this.nUpperBin; iFreq++) {
                    for (int iRow = 0; iRow < this.sensors; iRow++) {
                        for (int iCol = 0; iCol < this.sensors; iCol++) {
                            //this.psd[iRow][iCol][iFreq] = this.spec[iRow][iFreq]
                            //        .times(this.spec[iCol][iFreq].conjugate());

                            this.psd[iRow][iCol][iFreq] = new Complex(
                                this.spec[iRow][iFreq].re*this.spec[iCol][iFreq].re + this.spec[iRow][iFreq].im*this.spec[iCol][iFreq].im,
                                this.spec[iRow][iFreq].im*this.spec[iCol][iFreq].re - this.spec[iRow][iFreq].re*this.spec[iCol][iFreq].im);
                            //this.psd[iRow][iCol][iFreq].re = this.spec[iRow][iFreq].re*this.spec[iCol][iFreq].re + this.spec[iRow][iFreq].im*this.spec[iCol][iFreq].im;    
                            //this.psd[iRow][iCol][iFreq].im = this.spec[iRow][iFreq].im*this.spec[iCol][iFreq].re - this.spec[iRow][iFreq].re*this.spec[iCol][iFreq].im;
                            
                        }
                    }
                }
            } else {

                for (int iFreq = this.nLowerBin; iFreq < this.nUpperBin; iFreq++) {
                    for (int iRow = 0; iRow < this.sensors; iRow++) {
                        for (int iCol = 0; iCol < this.sensors; iCol++) {
                            //this.psd[iRow][iCol][iFreq] = oneMinusAlphaCpx
                            //        .times(this.spec[iRow][iFreq].times(this.spec[iCol][iFreq].conjugate()))
                            //        .plus(alphaCpx.times(this.psd[iRow][iCol][iFreq]));

                            this.psd[iRow][iCol][iFreq].re = this.oneMinusAlpha * 
                                (this.spec[iRow][iFreq].re*this.spec[iCol][iFreq].re + this.spec[iRow][iFreq].im*this.spec[iCol][iFreq].im) +
                                this.alpha * this.psd[iRow][iCol][iFreq].re;
                            this.psd[iRow][iCol][iFreq].im = this.oneMinusAlpha * 
                                (this.spec[iRow][iFreq].im*this.spec[iCol][iFreq].re - this.spec[iRow][iFreq].re*this.spec[iCol][iFreq].im) +
                                this.alpha * this.psd[iRow][iCol][iFreq].im;        
                        }
                    }
                }
            }


            // Testing of various azimuth angles

            //Complex[][] vectorA = new Complex[this.sensors][blocklen_half];
            double nP = 0;

            for (int iTheta = 0; iTheta < num_theta; iTheta++) {

                nP = 0;

                //double n_theta = theta[iTheta];

                // Calculation of Delays
                //double[] v_tau_theta = calculateDelays(n_theta, coordinates);


                Complex[][] vectorA = getSteeringVector(iTheta);

                // Creation of Steering Vector
                /*for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                    for (int iFreq = 0; iFreq < blocklen_half; iFreq++) {
                        vectorA[iSensor][iFreq] = (new Complex(0,
                                -2f * Math.PI * frequencies[iFreq] * v_tau_theta[iSensor])).exp();
                    }
                }*/

                

                // Power formula 27
                for (int iBin = this.nLowerBin; iBin < this.nUpperBin; iBin++) {

                    Complex trace = getRealTrace(this.psd, iBin);
                
                    Complex[][] eyeTimesTraceMinusPSD = generateEyeTimesComplex(this.sensors, trace);

                    for (int iRow = 0; iRow < this.sensors; iRow++) {
                        for (int iCol = 0; iCol < this.sensors; iCol++) {

                            eyeTimesTraceMinusPSD[iRow][iCol].re -= this.psd[iRow][iCol][iBin].re;
                            eyeTimesTraceMinusPSD[iRow][iCol].im -= this.psd[iRow][iCol][iBin].im;

                            //eyeTimesTraceMinusPSD[iRow][iCol] = eyeTimesTraceMinusPSD[iRow][iCol]
                            //        .minus(psd[iRow][iCol][iBin]);
                        }
                    }

                    //Complex[] tmp = this.complexSensorVector.clone();
                    Complex[] tmp = new Complex[this.sensors];
                    for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                        tmp[iSensor] = new Complex(0,0);
                    }

                    for (int iRow = 0; iRow < this.sensors; iRow++) {
                        for (int iCol = 0; iCol < this.sensors; iCol++) {
                            //tmp[iRow] = tmp[iRow].plus(vectorA[iCol][iBin].times(eyeTimesTraceMinusPSD[iRow][iCol]));

                            tmp[iRow].re += (eyeTimesTraceMinusPSD[iRow][iCol].re * vectorA[iCol][iBin].re - 
                                                eyeTimesTraceMinusPSD[iRow][iCol].im * vectorA[iCol][iBin].im); 
                            tmp[iRow].im += (eyeTimesTraceMinusPSD[iRow][iCol].im * vectorA[iCol][iBin].re + 
                                                eyeTimesTraceMinusPSD[iRow][iCol].re * vectorA[iCol][iBin].im); 
                        }
                    }

                    double tmp2 = 0;
                    for (int iRow = 0; iRow < this.sensors; iRow++) {
                        for (int iCol = 0; iCol < this.sensors; iCol++) {

                            //tmp2 += (vectorA[iCol][iRow].times(tmp[iRow])).abs();
                            tmp2 += Math.sqrt(Math.pow(vectorA[iCol][iBin].re*tmp[iRow].re - 
                                vectorA[iCol][iBin].im*tmp[iRow].im, 2) + 
                                Math.pow(vectorA[iCol][iBin].re*tmp[iRow].im + 
                                vectorA[iCol][iBin].im*tmp[iRow].re, 2)); 
// SQRT needed?
                        }
                    }

                    // Power summation over all (interesting) bins
                    nP += tmp2;
                }
                
                v_P_abs_sum[iTheta] = 1f / nP;
                mAbs_sum[iBlock][iTheta] = 1f / nP;

            }

            /*System.out.println("BLOICK:");

            for (int iTheta = 0; iTheta < v_P_abs_sum.length; iTheta++) {
                System.out.println(v_P_abs_sum[iTheta] + ",");
            }*/

            int arg = argMax(v_P_abs_sum);
            angles[iBlock] = theta[arg];
            
        }

    return angles;
}

}