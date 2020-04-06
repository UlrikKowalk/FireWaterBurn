
import java.io.File;


class HRTF{

    private String filename = "hrir/CH07IK25_0.mat";
    private double[][][] mHRTF;
    private MatReader matReader;

    public HRTF() {

        this.matReader = new MatReader(new File(this.filename));
        readHRTF();
    }

    private void readHRTF() {

        this.matReader.readAll();
        

         

    }

    public double[][] getHRTF(double direction) {

        double[][] tmp = new double[0][0];
        return tmp;
    }




}