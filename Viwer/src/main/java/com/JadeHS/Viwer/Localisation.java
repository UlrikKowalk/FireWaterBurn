package com.jadehs.viwer;

class Localisation {

    private int sensors;
    private int blocksize;
    private Complex[][][] psd;
    private Complex[][][] mTheta;
    private double alpha;
    private double oneMinusAlpha;
    private int nLowerBin, nUpperBin;
    private int blocksize_half;
    private int num_theta;
    private double[] v_P_abs_sum;
    private boolean firstBlock = true;
    private double[] theta;
    private SteeringVector steeringVector;
    private int nPeakImportance = 3;
    private boolean use_imaginary = true;

    public Localisation(int sensors, int samplerate, int blocksize, String arrayName) {

        this.sensors = sensors;
        this.blocksize = blocksize;

        //this.nLowerBin = 100;
        this.nLowerBin = 55;
        this.blocksize_half = (int) (this.blocksize / 2 + 1);
        this.nUpperBin = this.blocksize_half;
        
        this.alpha = 0.4f;
        this.oneMinusAlpha = 1.0f - this.alpha;

        double max_theta = 180;
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

        this.steeringVector = new SteeringVector(sensors, samplerate, blocksize, arrayName);
        this.mTheta = this.steeringVector.generateDelayTensor_DSB(this.theta);

        this.v_P_abs_sum = new double[this.num_theta];

        

        System.out.println("Theta: " + this.theta.length);

        // Power summation over all (interesting) bins
        if (this.use_imaginary) {
            System.out.println("Using imaginary.");
        } else {
            System.out.println("Only using real valued.");
        }

    }

    public int getNumTheta() {
        return this.num_theta;
    }

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
                        if (this.use_imaginary) {                    
                            tmp[iRow].im += (eyeTimesTraceMinusPSD[iRow][iCol].im * vectorA[iCol][iBin].re + 
                                                eyeTimesTraceMinusPSD[iRow][iCol].re * vectorA[iCol][iBin].im); 
                        }
                    }
                }

                Complex tmp2 = new Complex(0, 0);
                for (int iRow = 0; iRow < this.sensors; iRow++) {
                    for (int iCol = 0; iCol < this.sensors; iCol++) {

                        tmp2.re += vectorA[iCol][iBin].re*tmp[iRow].re + 
                            vectorA[iCol][iBin].im*tmp[iRow].im;
                        if (this.use_imaginary) {
                            tmp2.im += vectorA[iCol][iBin].re*tmp[iRow].im - 
                                vectorA[iCol][iBin].im*tmp[iRow].re; 
                        }

                    }
                }

                // Power summation over all (interesting) bins
                if (this.use_imaginary) {
                    nP += tmp2.abs();
                } else {
                    nP += tmp2.re;
                }
            }
            
            v_P_abs_sum[iTheta] = -1f / nP;

        }

        // Loudness Normalization
        /*double rms = Loudness.rms(v_P_abs_sum);
        for (int iTheta = 0; iTheta < this.num_theta; iTheta++) {
            v_P_abs_sum[iTheta] /= rms;
        }*/

        // Peak Picking
        int[] vPeaks = FindPeaks.findPeaks(v_P_abs_sum, nPeakImportance);

        // Quadratic Interpolation of Peak Positions
        double[][] mRealPeaks = QuadraticInterpolation.findRealPeaks(v_P_abs_sum, vPeaks);

        if (this.firstBlock) {
            this.firstBlock = false;
        }
        
        return mRealPeaks;
    }

}