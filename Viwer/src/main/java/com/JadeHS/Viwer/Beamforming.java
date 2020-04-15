package com.JadeHS.Viwer;

class Beamforming {

    private int sensors;
    private int samplerate;
    private int blocksize;
    private int fftlenHalf;
    private double[][] coordinates;
    private SteeringVector steeringvetor;
    private HRTF hrtf;

    public Beamforming(int sensors, int samplerate, int blocksize, double[][] coordinates, int num_theta_loc) {

        this.sensors = sensors;
        this.samplerate = samplerate;
        this.blocksize = blocksize;
        this.fftlenHalf = (int) blocksize / 2 + 1;
        this.coordinates = coordinates;

        this.steeringvetor = new SteeringVector(this.sensors, this.samplerate, this.blocksize);

        this.hrtf = new HRTF(num_theta_loc);


    }


    public double[][] filter(Complex[][] spec, double[] directions) {

        // Bins, L/R, Real/Imag
        double[][][] stereoSpec = new double[this.fftlenHalf][2][2];
        // Samples, L/R
        double[][] stereoOut = new double[this.blocksize][2];
      
        for (int iDir = 0; iDir < directions.length; iDir++) {

            // Sensors, Bins, Real/Imag
            double[][][] mTheta = this.steeringvetor.generateDelayTensor(directions[iDir]);
            // Bin, L/R , Real/Imag
            double[][][] hrtfCoefficients = this.hrtf.getHRTF(directions[iDir]);
            // Bins, Real/Imag
            double[][] monoSpec = new double[this.fftlenHalf][2];
            for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
                monoSpec[iBin][0] = 0.0f;
                monoSpec[iBin][1] = 0.0f;
            }

            // Beamforming
            for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
                    monoSpec[iBin][0] += spec[iSensor][iBin].re * mTheta[iSensor][iBin][0] - spec[iSensor][iBin].im * mTheta[iSensor][iBin][1];
                    monoSpec[iBin][1] += spec[iSensor][iBin].re * mTheta[iSensor][iBin][1] + spec[iSensor][iBin].im * mTheta[iSensor][iBin][0];
                }
            }

            

            // HRTFs
            for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
                // Real, Left
                stereoSpec[iBin][0][0] += monoSpec[iBin][0] * hrtfCoefficients[iBin][0][0] - monoSpec[iBin][1] * hrtfCoefficients[iBin][0][1];
                // Imag, Left
                stereoSpec[iBin][0][1] += monoSpec[iBin][0] * hrtfCoefficients[iBin][0][1] + monoSpec[iBin][1] * hrtfCoefficients[iBin][0][0];
                // Real, Right
                stereoSpec[iBin][1][0] += monoSpec[iBin][0] * hrtfCoefficients[iBin][1][0] - monoSpec[iBin][1] * hrtfCoefficients[iBin][1][1];
                // Imag, Right
                stereoSpec[iBin][1][1] += monoSpec[iBin][0] * hrtfCoefficients[iBin][1][1] + monoSpec[iBin][1] * hrtfCoefficients[iBin][1][0];
            }
        }

        // IFFT
        Complex[] tmp = new Complex[this.blocksize];

        // Left
        for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
            tmp[iBin] = new Complex(stereoSpec[iBin][0][0], stereoSpec[iBin][0][1]);
        }
        for (int iBin = 1; iBin < this.fftlenHalf-1; iBin++) {
            tmp[blocksize - iBin] = new Complex(stereoSpec[iBin][0][0], -stereoSpec[iBin][0][1]);
        }
        tmp = FFT.ifft(tmp);
        for (int iSample = 0; iSample < this.blocksize; iSample++) {
            stereoOut[iSample][0] = tmp[iSample].re;
        }

        // Right
        for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
            tmp[iBin] = new Complex(stereoSpec[iBin][1][0], stereoSpec[iBin][1][1]);
        }
        for (int iBin = 1; iBin < this.fftlenHalf-1; iBin++) {
            tmp[blocksize - iBin] = new Complex(stereoSpec[iBin][1][0], -stereoSpec[iBin][1][1]);
        }
        tmp = FFT.ifft(tmp);
        for (int iSample = 0; iSample < this.blocksize; iSample++) {
            stereoOut[iSample][1] = tmp[iSample].re;
        }

        return stereoOut;
    }

    


}