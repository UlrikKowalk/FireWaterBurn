package com.jadehs.viwer;

class QuadraticInterpolation {

    static double[] quadraticInterpolation(double magPrev, double mag, double magNext) {

        double[] result = new double[2];
        result[0] = (magNext - magPrev) / (2.0 * (2.0 * mag - magNext - magPrev));
        result[1] = mag - 0.25 * (magPrev - magNext) * result[0];
        return result;

    }     
    
    static double[][] findRealPeaks(double[] inArray, int[] vPeaks) {

        double[][] vRealPeaks = new double[vPeaks.length][2];

        for (int iPeak = 0; iPeak < vPeaks.length; iPeak++) {

            double[] tmp = new double[2];

            // Avoid circular confusion
            if (vPeaks[iPeak] == 0) {
                tmp = quadraticInterpolation(inArray[inArray.length-1], inArray[0], inArray[1]);
            } else if (vPeaks[iPeak] == inArray.length-1) {
                tmp = quadraticInterpolation(inArray[inArray.length-2], inArray[inArray.length-1], inArray[0]);
            } else {
                tmp = quadraticInterpolation(inArray[vPeaks[iPeak]-1], inArray[vPeaks[iPeak]], inArray[vPeaks[iPeak]+1]);
            }

            vRealPeaks[iPeak][0] = tmp[0] + vPeaks[iPeak];
            vRealPeaks[iPeak][1] = tmp[1];
        }

        return vRealPeaks;
    }

}