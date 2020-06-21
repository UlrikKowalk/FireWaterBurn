package com.jadehs.viwer;

import java.io.File;
import com.jadehs.viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.io.*;
import com.jadehs.viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.types.*;


class HRTF {

    private String filename_hrtf = "hrir/CH07IK25_0.mat";
    private String filename_directions = "hrir/directions.mat";
    private double[][][] mHRTF_L;
    private double[][][] mHRTF_R;
    private MatFileReader matReader;
    private File file;
    private int num_theta;
    private int num_theta_loc;
    private int fftlenHalf;
    private int nChannels = 2;
    private double[] directions;
    private double theta;
    private double thetaTransform;
    private double thetaFactor;
    

    public HRTF(int num_theta_loc) {

        // Obtain measured hrtf directions
        this.file = new File(this.filename_directions);
        readDirections();

        // Obtain measured hrtfs
        this.file = new File(this.filename_hrtf);
        readHRTF();

        this.num_theta_loc = num_theta_loc;
        this.thetaTransform = this.theta * this.num_theta / this.num_theta_loc;
        this.thetaFactor = this.num_theta / this.num_theta_loc;

    }

    private void readHRTF() {

        try {

            this.matReader = new MatFileReader(this.file);

            MLDouble data = (MLDouble) this.matReader.getMLArray("mHRTF_L");
            this.fftlenHalf = data.getN();

            this.mHRTF_L = new double[this.num_theta][this.fftlenHalf][2];
            for (int iDir = 0; iDir < this.num_theta; iDir++) {
                for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
                    this.mHRTF_L[iDir][iBin][0] = data.getReal(iDir, iBin);
                    this.mHRTF_L[iDir][iBin][1] = data.getImaginary(iDir, iBin);
                }
            }

            this.mHRTF_R = new double[this.num_theta][this.fftlenHalf][2];
            data = (MLDouble) this.matReader.getMLArray("mHRTF_R");
            for (int iDir = 0; iDir < this.num_theta; iDir++) {
                for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
                    this.mHRTF_R[iDir][iBin][0] = data.getReal(iDir, iBin);
                    this.mHRTF_R[iDir][iBin][1] = data.getImaginary(iDir, iBin);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }

    private void readDirections() {

        try {

            this.matReader = new MatFileReader(this.file);
            MLDouble data = (MLDouble) this.matReader.getMLArray("vDirections");
            this.num_theta = data.getM();
            this.directions = new double[this.num_theta];
            for (int iDir = 0; iDir < this.num_theta; iDir++) {
                this.directions[iDir] = data.get(iDir);
            }
            this.theta = 360.f / this.num_theta;

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public double[][][] getHRTF(double direction) {

        // Bin, L/R, Real/Imag
        double[][][] tmp = new double[this.fftlenHalf][2][this.nChannels];

        int[] directionIdx = getDirectionIndex(direction);
        double factorOfFirst = getDirectionFactor(direction, directionIdx);

        for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
            // Real, Left
            tmp[iBin][0][0] = factorOfFirst * this.mHRTF_L[directionIdx[0]][iBin][0] + (1.0f - factorOfFirst) * this.mHRTF_L[directionIdx[1]][iBin][0];
            // Imag, Left
            tmp[iBin][0][1] = factorOfFirst * this.mHRTF_L[directionIdx[0]][iBin][1] + (1.0f - factorOfFirst) * this.mHRTF_L[directionIdx[1]][iBin][1];
            // Real, Right
            tmp[iBin][1][0] = factorOfFirst * this.mHRTF_R[directionIdx[0]][iBin][0] + (1.0f - factorOfFirst) * this.mHRTF_R[directionIdx[1]][iBin][0];
            // Imag, Right
            tmp[iBin][1][1] = factorOfFirst * this.mHRTF_R[directionIdx[0]][iBin][1] + (1.0f - factorOfFirst) * this.mHRTF_R[directionIdx[1]][iBin][1];
            
            /*
            // Real, Left
            tmp[iBin][0][0] = this.mHRTF_L[directionIdx[0]][iBin][0];
            // Imag, Left
            tmp[iBin][0][1] = this.mHRTF_L[directionIdx[0]][iBin][1];
            // Real, Right
            tmp[iBin][1][0] = this.mHRTF_R[directionIdx[0]][iBin][0];
            // Imag, Right
            tmp[iBin][1][1] = this.mHRTF_R[directionIdx[0]][iBin][1];
            */
        }

        return tmp;
    }

    private int[] getDirectionIndex(double angle) {
        
        int minIdx = (int) ( angle * this.thetaFactor);

        int[] tmp = new int[2];
        if (minIdx == this.num_theta-1) {
            // Cyclic value domain
            tmp = new int[] {this.num_theta-1, 0};
        } else {
            tmp = new int[] {minIdx, minIdx+1};
        }

        return tmp;
    }

    private double getDirectionFactor(double angle, int[] idx) {

        return 1.0f - (angle * this.thetaTransform - (double) directions[idx[0]]) / this.theta;

    }

}