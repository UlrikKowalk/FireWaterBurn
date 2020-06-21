package com.jadehs.viwer;

class Loudness {

    static double rms(double[] inArray) {
        
        double len = inArray.length;
        double sum = 0.0f;

        for (int iSample = 0; iSample < len; iSample++) {
            sum += inArray[iSample] * inArray[iSample];
        }
        sum /= len;

        double rms = Math.sqrt(sum);

        return rms;
    }

    static void rmsTest() {
        double[] testArray = new double[] {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f};
        double expectedResult = 6.2048f;

        if (Math.abs(rms(testArray) - expectedResult) < 0.001) {
            System.out.println("Test successful.");
        } else {
            System.out.println("Test not successful.");
        }

    }

    public static void main(String[] args) {
        rmsTest();
    }

}