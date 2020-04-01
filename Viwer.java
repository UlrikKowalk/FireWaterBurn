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
    private Complex[][][] psd;

    private Localisation localisation;
    private SourceManager sourceManager;
    private SteeringVector steeringVector;
    private Beamforming beamforming;

    public Viwer(int samplerate, int blocksize, double overlap, int sensors) {

        this.samplerate = samplerate;
        this.blocksize = blocksize;
        this.overlap = overlap;
        this.sensors = sensors;

        this.hopsize = (int) (this.blocksize * (1.0 - overlap));
        this.blocklen_half = (int) (this.blocksize / 2 + 1);
        
        this.steeringVector = new SteeringVector(this.sensors, this.samplerate, blocksize);

        this.localisation = new Localisation(this.sensors, this.samplerate, this.blocksize);
        this.num_theta = this.localisation.getNumTheta();

        this.dt = Kalman.calculateDt(this.blocksize, this.samplerate);

        this.sourceManager = new SourceManager(dt, this.num_theta);

        // Pre-Allocation and Initialisation of signal arrays
        this.signalBlock = new Complex[this.blocksize];
        for (int iSample = 0; iSample < this.blocksize; iSample++) {
            this.signalBlock[iSample] = new Complex(0,0);
        }
        this.spec = new Complex[this.sensors][this.blocklen_half];
        this.psd = new Complex[this.sensors][this.sensors][this.blocklen_half];

    }

    public double[][] analyseAndFilter(double[][] sensor_data) {

        int length = sensor_data[0].length;
        int num_blocks = (int) Math.floor((length - blocksize) / hopsize);

        double[] window = hann(this.blocksize);

        double[][] mAccu = new double[num_blocks][10];
        
        int idx_in, idx_out;
        
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
                }
            }

            // Localisation
            double[][] directions = this.localisation.sourceDirections(this.spec);

            // Tracking
            double[] vKalmanPeaks = this.sourceManager.trackSources(directions);




            
            // Write results to matrix
            for (int iSource = 0; iSource < vKalmanPeaks.length; iSource++) {
                mAccu[iBlock][iSource] = vKalmanPeaks[iSource];
            }
            for (int iSource = vKalmanPeaks.length; iSource < SourceManager.MAX_SOURCES; iSource++) {
                mAccu[iBlock][iSource] = -255.0f;
            }
     
        }







        // Write results to file
        ResultWriter resultWriter = new ResultWriter("threshold.txt");
        for (int iBlock = 0; iBlock < num_blocks; iBlock++) {
            resultWriter.write(mAccu[iBlock]);
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