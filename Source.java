
public class Source {

    private SourceManager context;
    private Kalman kalman;
    private double nAngle;
    private double nSpread;
    private double nWidth = 4.0f;
    private int blocksWithoutUpdate = 0;
    private double alpha = 0.95f;

    public Source(double angle, SourceManager sourceManager, double dt, double nSpread) {

        this.nAngle = angle;
        this.context = sourceManager;
        this.nSpread = nSpread;
        this.kalman = new Kalman(dt, this.nAngle);

    }


    public double[] iterate(double nCandidate) {
        this.blocksWithoutUpdate = 0;
        this.nAngle = (1.0f - this.alpha) * this.kalman.iterate(nCandidate) + this.alpha * this.nAngle;
        return new double[] {this.nAngle, Math.pow(0.5f, this.blocksWithoutUpdate)};
    }

    public double[] iterateWeighted(double[] vCandidates) {

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
            tmp = pdf(vTempCandidates[iWeight], this.nAngle, this.nWidth);
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
        return kalman.getEstimate();
    }

    public int getBlocksWithoutUpdate() {
        return this.blocksWithoutUpdate;
    }

    public double getAngle() {
        return this.kalman.getEstimate();
    }

    public double getPdf(double x) {
        return pdf(x, getAngle(), this.nWidth);
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