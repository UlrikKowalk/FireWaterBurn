package com.JadeHS.Viwer;

class Beamforming {

    private int sensors;
    private int samplerate;
    private int blocksize;
    private double[][] coordinates;
    private SteeringVector steeringvetor;
    private HRTF hrtf;

    public Beamforming(int sensors, int samplerate, int blocksize, double[][] coordinates) {

        this.sensors = sensors;
        this.samplerate = samplerate;
        this.blocksize = blocksize;
        this.coordinates = coordinates;

        this.steeringvetor = new SteeringVector(this.sensors, this.samplerate, this.blocksize);

        this.hrtf = new HRTF();


    }


    private double[][] filter(double[][] spec, double[] directions) {

        // TODO: Include HRTFs at this step
      
        for (int iDir = 0; iDir < directions.length; iDir++) {

            double[][] hrtfCoefficients = this.hrtf.getHRTF(directions[iDir]);
            //double[][][] mTheta = SteeringVector.generate 


        }





        double[][] tmp = new double[0][0];
        return tmp;
    }

    


}