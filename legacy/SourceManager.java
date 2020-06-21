import java.util.ArrayList;

public class SourceManager {

    public final static int MAX_SOURCES = 5;
    public final static int MAX_BLOCKS_WITHOUT_UPDATE = 5;
    private int NUM_SOURCES = 0;

    private ArrayList<Source> sources;
    private ArrayList<ArrayList<Double>> candidateList;

    private double dt;
    private double num_theta;
    private double thresholdProbability = 0.05f;
    private double thresholdPeak = 1.8f;

    private double[] vKalmanPeaks;
    private double[] vKalmanPeaksWeighted;
    private double[] vCandidates;

    public SourceManager(double dt, double num_theta) {

        this.dt = dt;
        this.num_theta = num_theta;
        this.sources = new ArrayList<>();
        
        this.vKalmanPeaks = new double[SourceManager.MAX_SOURCES];
        this.vKalmanPeaksWeighted = new double[SourceManager.MAX_SOURCES];
        for (int iSource = 0; iSource < SourceManager.MAX_SOURCES; iSource++) {
            this.vKalmanPeaks[iSource] = -255;
            this.vKalmanPeaksWeighted[iSource] = -255;
        }
    }

    public double[] trackSources(double[][] mPeaks) {

        this.candidateList = new ArrayList<>();
        ArrayList<Integer> sourcesToDelete = new ArrayList<>();

        for (int iSource = 0; iSource < this.NUM_SOURCES; iSource++) {
            this.candidateList.add(new ArrayList<>());
        }

        for (int iPeak = 0; iPeak < mPeaks[0].length; iPeak++) {

            if (this.NUM_SOURCES == 0) {
                // No sources yet -> every peak is a potential new source
                double[] firstResult = new double[2 * SourceManager.MAX_SOURCES];
                int idxtmp = 0;
                for (int iPk = 0; iPk < mPeaks[0].length; iPk++) {
                    if (mPeaks[iPk][1] > this.thresholdPeak && this.NUM_SOURCES < SourceManager.MAX_SOURCES) {
                        addSource(mPeaks[iPk][0]);
                        firstResult[idxtmp] = mPeaks[iPk][0];
                        double[] tmp = this.sources.get(this.sources.size()-1).iterate(mPeaks[iPk][0]);
                        this.vKalmanPeaks[idxtmp] = tmp[0];
                        idxtmp++;
                    }
                }
                return firstResult;

            } else {

                if (mPeaks[iPeak][1] > this.thresholdPeak) {

                    // Sources exist
                    double[] vProbabilities = getProbabilities(mPeaks[iPeak][0]);
                    int idx = argMax(vProbabilities, this.thresholdProbability);

                    if (idx == -255) {
                        // No matching source found
                        if (this.NUM_SOURCES < SourceManager.MAX_SOURCES) {
                            addSource(mPeaks[iPeak][0]);
                            this.candidateList.get(this.candidateList.size()-1).add(mPeaks[iPeak][0]); 
                        }
                    } else {
                        // Matching source found
                        this.candidateList.get(idx).add(mPeaks[iPeak][0]); 
                    }
                }
            }
        }

        // Kalman Filtering
        for (int iSource = 0; iSource < this.NUM_SOURCES; iSource++) {
            
            if (candidateList.size() > 0) {
                this.vCandidates = new double[candidateList.get(iSource).size()];
                for (int iCandidate = 0; iCandidate < candidateList.get(iSource).size(); iCandidate++) {
                    this.vCandidates[iCandidate] = this.candidateList.get(iSource).get(iCandidate);
                }
            } else {
                this.vCandidates = new double[] {mPeaks[iSource][0], mPeaks[iSource][0]};
            }
          
            if (candidateList.get(iSource).size() == 0) {
                this.vKalmanPeaks[iSource] = this.sources.get(iSource).noUpdate();
            } else if (candidateList.get(iSource).size() == 1) {
                double[] tmp = this.sources.get(iSource).iterate(vCandidates[0]);
                this.vKalmanPeaks[iSource] = tmp[0];
            } else {
                double[] tmp = this.sources.get(iSource).iterateWeighted(vCandidates);
                this.vKalmanPeaks[iSource] = tmp[0];
            }
             
            // Sources with more blocks without update than MAX_BLOCKS_WITHOUT_UPDATE are deleted
            if (this.sources.get(iSource).getBlocksWithoutUpdate() >= SourceManager.MAX_BLOCKS_WITHOUT_UPDATE) {
                sourcesToDelete.add(iSource);
            }
        }

        deleteSource(sourcesToDelete);

        double[] vResult = new double[2 * SourceManager.MAX_SOURCES];
        for (int iSource = 0; iSource < this.NUM_SOURCES; iSource++) {
            vResult[iSource] = this.vKalmanPeaks[iSource];
        }
        for (int iSource = this.NUM_SOURCES; iSource < 2 * SourceManager.MAX_SOURCES; iSource++) {
            vResult[iSource] = -255.0f;
        }
        return vResult;
    }

    private void addSource(double nPeak) {
        this.sources.add(new Source(nPeak, this, this.dt, this.num_theta));
        this.NUM_SOURCES++;
        this.candidateList.add(new ArrayList<>());
    }

    private void deleteSource(ArrayList<Integer> delete) {
        if (delete.size() > 0) {
            for (int iSource = delete.size()-1; iSource >= 0; iSource--) {
                this.sources.remove((int) delete.get(iSource));
                this.NUM_SOURCES--;
            }
        }
    }

    private double[] getProbabilities(double nCandidate) {
        double[] vProbabilities = new double[this.sources.size()]; 
        for (int iSource = 0; iSource < this.sources.size(); iSource++) {
            vProbabilities[iSource] = this.sources.get(iSource).getPdf(nCandidate);
        }
        return vProbabilities;
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

    private int argMax(double[] inVec, double threshold) {
        int arg = -255;
        double max = 0;
        for (int idx = 0; idx < inVec.length; idx++) {
            if (inVec[idx] > max && inVec[idx] > threshold) {
                max = inVec[idx];
                arg = idx;
            }
        }
        return arg;
    }
}