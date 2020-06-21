package com.jadehs.viwer;

public class Array {

    private String name;
    private double[][] coordinates;
    private int sensors;

    public Array(String name) {
        this.name = name;

        switch(this.name) {
            case "MatrixVoice":
                this.sensors = 8;
                establishCoordinates_MatrixVoice();
                break;
            case "Boomerang":
                this.sensors = 15;
                establishCoordinates_Boomerang();
                break;
        }
    }

    public String getName() {
        return this.name;
    }

    public double[][] getCoordinates() {
        return this.coordinates;
    }

    private void establishCoordinates_MatrixVoice() {

        //Generate Coordinates for MatrixVoice Array
        
        this.coordinates = new double[this.sensors][3];
        double omega;

        // Not necessary, just feels right
        this.coordinates[0][0] = 0f;
        this.coordinates[0][1] = 0f;
        this.coordinates[0][2] = 0f;

        for (int iSensor = 1; iSensor < this.sensors; iSensor++) {
            omega = -2f * Math.PI * (iSensor - 1) / (this.sensors - 1) - Math.PI;
            this.coordinates[iSensor][0] = Math.PI * Math.cos(omega) * 0.01f;
            this.coordinates[iSensor][1] = Math.PI * Math.sin(omega) * 0.01f;
            this.coordinates[iSensor][2] = 0f;
        }
    }

    private void establishCoordinates_Boomerang() {

        // Generate Coordinates for the Viwer Boomerang v.01

        this.coordinates = new double[this.sensors][3];
        this.coordinates[0] = new double[] {-0.11381f,-0.17055f,0f};
        this.coordinates[1] = new double[] {-0.07911f,-0.13218f,0f};
        this.coordinates[2] = new double[] {-0.052397f,-0.10086f,0f};
        this.coordinates[3] = new double[] {-0.03225f,-0.074911f,0f};
        this.coordinates[4] = new double[] {-0.017594f,-0.052949f,0f};
        this.coordinates[5] = new double[] {-0.0076496f,-0.033807f,0f};
        this.coordinates[6] = new double[] {-0.0018873f,-0.016465f,0f};
        this.coordinates[7] = new double[] {0f,0f,0f};
        this.coordinates[8] = new double[] {-0.0018873f,0.016465f,0f};
        this.coordinates[9] = new double[] {-0.0076496f,0.033807f,0f};
        this.coordinates[10] = new double[] {-0.017594f,0.052949f,0f};
        this.coordinates[11] = new double[] {-0.03225f,0.074911f,0f};
        this.coordinates[12] = new double[] {-0.052397f,0.10086f,0f};
        this.coordinates[13] = new double[] {-0.07911f,0.13218f,0f};
        this.coordinates[14] = new double[] {-0.11381f,0.17055f,0f};

    }

}