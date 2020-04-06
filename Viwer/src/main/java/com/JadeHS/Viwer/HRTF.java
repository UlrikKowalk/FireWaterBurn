package com.JadeHS.Viwer;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;

//import com.JadeHS.Viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.io.MatFileReader;
//import com.JadeHS.Viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.io.MatFileHeader;
import com.JadeHS.Viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.io.*;
import com.JadeHS.Viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.common.*;
import com.JadeHS.Viwer.external.JMatIO.src.main.java.ca.mjdsystems.jmatio.types.*;


class HRTF{

    private String filename_hrtf = "hrir/CH07IK25_0.mat";
    private String filename_directions = "hrir/directions.mat";
    private double[][][] mHRTF_L;
    private double[][][] mHRTF_R;
    private MatFileReader matReader;
    private File file;
    private int num_theta;
    private int fftlenHalf;
    private int nChannels = 2;
    private double[] directions;
    

    public HRTF() {

        // Obtain measured hrtf directions
        this.file = new File(this.filename_directions);
        readDirections();

        // Obtain measured hrtfs
        this.file = new File(this.filename_hrtf);
        readHRTF();

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

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public double[][][] getHRTF(double direction) {

        // Bin, Real/Imag, L/R 
        double[][][] tmp = new double[this.fftlenHalf][2][this.nChannels];

        int[] directionIdx = new int[2];
        double factorOfFirst = 1.0f;

        for (int iBin = 0; iBin < this.fftlenHalf; iBin++) {
            // Real, Left
            tmp[iBin][0][0] = factorOfFirst * this.mHRTF_L[iBin][directionIdx[0]][0] + (1.0f - factorOfFirst) * this.mHRTF_L[iBin][directionIdx[1]][0];
            // Imag, Left
            tmp[iBin][1][0] = factorOfFirst * this.mHRTF_L[iBin][directionIdx[0]][1] + (1.0f - factorOfFirst) * this.mHRTF_L[iBin][directionIdx[1]][1];
            // Real, Right
            tmp[iBin][0][1] = factorOfFirst * this.mHRTF_R[iBin][directionIdx[0]][0] + (1.0f - factorOfFirst) * this.mHRTF_R[iBin][directionIdx[1]][0];
            // Imag, Right
            tmp[iBin][1][1] = factorOfFirst * this.mHRTF_R[iBin][directionIdx[0]][1] + (1.0f - factorOfFirst) * this.mHRTF_R[iBin][directionIdx[1]][1];
        }

        return tmp;
    }

    public int[] getDirectionIndex(double angle) {




        int[] tmp = new int[2];
        return tmp;
    }

    public double getDirectionFactor(double angle) {



        return 1.0f;
    }



}