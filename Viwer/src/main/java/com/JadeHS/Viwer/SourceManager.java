package com.jadehs.viwer;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

class SourceManager {

    public final static int MAX_SOURCES = 5;
    public final static double MAX_SECONDS_WITHOUT_UPDATE = 2.0f;
    public int MAX_BLOCKS_WITHOUT_UPDATE;
    private int NUM_SOURCES = 0;

    private ArrayList<Source> sources;
    private CandidateList candidateList;

    // MatrixVoice: 2.0f, 0.05f, 1.8f


    private double dt;
    private double num_theta;
    private double thresholdProbability = 0.01f;//0.05f;
    private double thresholdPeak = -60.0f;//1.8f; // imag: 2.0f;

    private double[] vKalmanPeaks;
    private double[] vKalmanPeaksWeighted;
    private double[] vCandidates;

    public SourceManager(double dt, double num_theta) {

        this.dt = dt;
        this.num_theta = num_theta;
        this.sources = new ArrayList<>();
        this.MAX_BLOCKS_WITHOUT_UPDATE = (int) (SourceManager.MAX_SECONDS_WITHOUT_UPDATE / this.dt);
        
        this.vKalmanPeaks = new double[SourceManager.MAX_SOURCES];
        this.vKalmanPeaksWeighted = new double[SourceManager.MAX_SOURCES];
        for (int iSource = 0; iSource < SourceManager.MAX_SOURCES; iSource++) {
            this.vKalmanPeaks[iSource] = -255;
            this.vKalmanPeaksWeighted[iSource] = -255;
        }
    }

    public double[] trackSources(double[][] mPeaks) {

        mPeaks = sortPeaks(mPeaks);

        this.candidateList = new CandidateList();
        ArrayList<Integer> sourcesToDelete = new ArrayList<>();

        int nPeaks = Math.min(SourceManager.MAX_SOURCES, mPeaks.length);

        for (int iSource = 0; iSource < this.NUM_SOURCES; iSource++) {
            this.candidateList.add(new ArrayList<Double>());
        }

        //for (int iPeak = 0; iPeak < mPeaks.length; iPeak++) {
        for (int iPeak = 0; iPeak < nPeaks; iPeak++) {

            if (this.NUM_SOURCES == 0) {
                // No sources yet -> every peak is a potential new source
                double[] tmpResult = new double[Math.min(mPeaks.length, SourceManager.MAX_SOURCES)];
                int idxtmp = 0;
                for (int iPk = 0; iPk < mPeaks.length; iPk++) {
                    if (mPeaks[iPk][1] > this.thresholdPeak && this.NUM_SOURCES < SourceManager.MAX_SOURCES) {
                        addSource(mPeaks[iPk][0]);
                        tmpResult[idxtmp] = mPeaks[iPk][0];
                        double tmp = this.sources.get(this.sources.size()-1).iterate(mPeaks[iPk][0]);
                        this.vKalmanPeaks[idxtmp] = tmp;
                        idxtmp++;
                    }
                }

                double[] firstResult = new double[idxtmp];
                for (int iPk = 0; iPk < idxtmp; iPk++){
                    firstResult[iPk] = tmpResult[iPk];
                }
                firstResult = new double[0]; // Fist time a source is seen -> no notice, let's see if it holds

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
                // Remove spurious peaks
                if (this.sources.get(iSource).getIsNew()) {
                    sourcesToDelete.add(iSource);
                } else {
                    this.vKalmanPeaks[iSource] = this.sources.get(iSource).noUpdate();
                }
            } else if (candidateList.get(iSource).size() == 1) {
                double tmp = this.sources.get(iSource).iterate(vCandidates[0]);
                this.vKalmanPeaks[iSource] = tmp;
            } else {
                double tmp = this.sources.get(iSource).iterateWeighted(vCandidates);
                this.vKalmanPeaks[iSource] = tmp;
            }
             
            // Sources with more blocks without update than MAX_BLOCKS_WITHOUT_UPDATE are deleted
            if (this.sources.get(iSource).getBlocksWithoutUpdate() >= this.MAX_BLOCKS_WITHOUT_UPDATE) {
                sourcesToDelete.add(iSource);
            }
        }

        deleteSource(sourcesToDelete);

        double[] vResult = new double[this.NUM_SOURCES];
        for (int iSource = 0; iSource < this.NUM_SOURCES; iSource++) {
            vResult[iSource] = this.vKalmanPeaks[iSource];
        }

        return vResult;
    }

    private double[][] sortPeaks(double[][] mPeaks) {
        java.util.Arrays.sort(mPeaks, new java.util.Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                return Double.compare(a[1], b[1]);
            }
        });
        return mPeaks;
    }

    public int getNumSources() {
        return this.NUM_SOURCES;
    }

    private void addSource(double nPeak) {
        this.sources.add(new Source(nPeak, this.dt, this.num_theta));
        this.NUM_SOURCES++;
        this.candidateList.add(new ArrayList<Double>());
    }

    private void deleteSource(ArrayList<Integer> delete) {
        if (delete.size() > 0) {
            for (int iSource = delete.size()-1; iSource >= 0; iSource--) {
                this.sources.remove((int) delete.get(iSource));
                this.NUM_SOURCES--;
            }
        }
    }

    private void deleteSource(Set<Integer> delete) {
        this.sources.removeAll(delete);
        this.NUM_SOURCES = this.sources.size();
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