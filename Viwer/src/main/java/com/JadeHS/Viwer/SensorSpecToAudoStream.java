package com.jadehs.viwer;

public class SensorSpecToAudoStream {


    public static double[] transformToMono(Complex[][] spec) {

        int blocksizehalf = spec[0].length;
        int sensors = spec.length;
        int blocksize = (blocksizehalf - 1) * 2;
        
        Complex[] tmp = new Complex[blocksize];
        for (int iBin = 0; iBin < blocksize; iBin++) {
            tmp[iBin] = new Complex(0.0f, 0.0f);
        }

        // Mono
        /*for (int iSensor = 0; iSensor < sensors; iSensor++) { 
            for (int iBin = 0; iBin < blocksizehalf; iBin++) {
                tmp[iBin] = tmp[iBin].plus(spec[iSensor][iBin]);
            }
            for (int iBin = 0; iBin < blocksizehalf-2; iBin++) {
                tmp[blocksize - iBin - 1] = tmp[blocksize - iBin - 1].plus(spec[iSensor][iBin]);
            }
        }*/

        for (int iSensor = 0; iSensor < sensors; iSensor++) { 
            for (int iBin = 0; iBin < blocksizehalf; iBin++) {
                tmp[iBin] = tmp[iBin].plus(spec[iSensor][iBin]);
            }
            for (int iBin = 1; iBin < blocksizehalf-1; iBin++) {
                tmp[blocksize - iBin] = tmp[blocksize - iBin].plus(spec[iSensor][iBin].conjugate());
            }
        }
    

        tmp = FFT.ifft(tmp);

        double[] monoOut = new double[blocksize];

        for (int iSample = 0; iSample < blocksize; iSample++) {
            monoOut[iSample] = tmp[iSample].re;
        }

        return monoOut;
    }


    public static double[] transformMono(Complex[] spec) {

        int blocksizehalf = spec.length;
        int blocksize = (blocksizehalf - 1) * 2;

        Complex[] tmp = new Complex[blocksize];

        for (int iBin = 0; iBin < blocksizehalf; iBin++) {
            tmp[iBin] = spec[iBin];
        }
        for (int iBin = 1; iBin < blocksizehalf-1; iBin++) {
            tmp[blocksize - iBin] = spec[iBin].conjugate();
        }

        tmp = FFT.ifft(tmp);

        double[] res = new double[blocksize];

        for (int iSample = 0; iSample < blocksize; iSample++) {
            res[iSample] = tmp[iSample].re;
        }

        return res;
    
    }
}