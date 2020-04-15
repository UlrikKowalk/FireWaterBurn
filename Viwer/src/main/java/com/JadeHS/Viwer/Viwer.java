package com.JadeHS.Viwer;

class Viwer {

    private int samplerate;
    private int blocksize;
    private int blocklen_half;
    private int hopsize;
    private int sensors;
    private int num_theta;
    private double overlap;
    private double dt;

    private Complex[] signalBlock;
    private Complex[][] spec;
    //private Complex[][][] psd;

    private Localisation localisation;
    private SourceManager sourceManager;
    private SteeringVector steeringVector;
    private Beamforming beamforming;

    public Viwer(int samplerate, int blocksize, double overlap, int sensors) {

        this.samplerate = samplerate;
        this.blocksize = blocksize;
        this.overlap = overlap;
        this.sensors = sensors;

        this.hopsize = (int) (this.blocksize * (1.0 - this.overlap));
        this.blocklen_half = (int) (this.blocksize / 2 + 1);
        
        this.steeringVector = new SteeringVector(this.sensors, this.samplerate, blocksize);

        this.localisation = new Localisation(this.sensors, this.samplerate, this.blocksize);
        this.num_theta = this.localisation.getNumTheta();

        this.dt = Kalman.calculateDt(this.hopsize, this.samplerate);

        this.sourceManager = new SourceManager(this.dt, this.num_theta);

        this.beamforming = new Beamforming(this.sensors, this.samplerate, this.blocksize, this.steeringVector.getCoordinates(), this.num_theta);

        // Pre-Allocation and Initialisation of signal arrays
        this.signalBlock = new Complex[this.blocksize];
        for (int iSample = 0; iSample < this.blocksize; iSample++) {
            this.signalBlock[iSample] = new Complex(0,0);
        }
        this.spec = new Complex[this.sensors][this.blocklen_half];
    }

    public double[][] analyseAndFilter(double[][] sensor_data) {

        int length = sensor_data[0].length;
        int num_blocks = (int) Math.floor((length - blocksize) / this.hopsize);

        double[] window = hann(this.blocksize);
        for (int iSample = 0; iSample < this.blocksize; iSample++) {
            window[iSample] = Math.sqrt(window[iSample]);
        }

        double[][] mAccu = new double[num_blocks][10];

        double[][] mOut = new double[length][2];
        
        int idx_in, idx_out = 0;
        
        for (int iBlock = 0; iBlock < num_blocks; iBlock++) {

            idx_in = (int) (iBlock * this.hopsize);
            idx_out = (int) (idx_in + blocksize);

            // Generation of the Spectrum
            for (int iSensor = 0; iSensor < this.sensors; iSensor++) {
                for (int iSample = idx_in; iSample < idx_out; iSample++) {
                    signalBlock[iSample - idx_in].re = sensor_data[iSensor][iSample] * window[iSample - idx_in];
                    signalBlock[iSample - idx_in].im = 0.0f;
                }
                InplaceFFT.fft(signalBlock);
                for (int iFreq = 0; iFreq < this.blocklen_half; iFreq++) {
                    this.spec[iSensor][iFreq] = new Complex(signalBlock[iFreq].re, signalBlock[iFreq].im);
                }
            }

            // Localisation
            double[][] directions = this.localisation.sourceDirections(this.spec);

            // Tracking
            double[] vKalmanPeaks = this.sourceManager.trackSources(directions);

            // Beamforming and Positioning
            double[][] mAudioOut = beamforming.filter(this.spec, vKalmanPeaks);

            for (int iSample = 0; iSample < this.blocksize; iSample++) {
                mOut[idx_in + iSample][0] = mOut[idx_in + iSample][0] + mAudioOut[iSample][0] * window[iSample];
                mOut[idx_in + iSample][1] = mOut[idx_in + iSample][1] + mAudioOut[iSample][1] * window[iSample];
            }


            // Write directions to matrix
            for (int iSource = 0; iSource < directions.length; iSource++) {
                mAccu[iBlock][iSource] = directions[iSource][0];
            }
            // Write Kalman filtered source directions to matrix
            for (int iSource = 0; iSource < vKalmanPeaks.length; iSource++) {
                mAccu[iBlock][SourceManager.MAX_SOURCES + iSource] = vKalmanPeaks[iSource];
            }
        
        }

        
        // Write directions (normal, filtered) to file
        ResultWriter resultWriter = new ResultWriter("threshold.txt");
        resultWriter.write(SourceManager.MAX_SOURCES);
        resultWriter.write(this.num_theta);
        for (int iBlock = 0; iBlock < num_blocks; iBlock++) {
            resultWriter.write(mAccu[iBlock]);
        }

        // Write audio data to file
        resultWriter = new ResultWriter("Viwer_Out.txt");
        for (int iSample = 0; iSample < idx_out; iSample++) {
            resultWriter.write(mOut[iSample][0]);
            resultWriter.write(mOut[iSample][1]);
        }

        double[][] tmp = new double[0][0];

        return tmp;
    }

    private double[] hann(int blocksize) {
        double[] window = new double[blocksize];
        for (int iSample = 0; iSample < blocksize; iSample++) {
            window[iSample] = 0.5 * (1 - Math.cos(2 * Math.PI * iSample / (blocksize - 1)));
        }
        return window;
    }


}