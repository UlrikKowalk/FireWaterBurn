
class QuadraticInterpolation {

    static double[] quadraticInterpolation(double magPrev, double mag, double magNext) {

        double[] result = new double[2];
        result[0] = (magNext - magPrev) / (2.0 * (2.0 * mag - magNext - magPrev));
        result[1] = mag - 0.25 * (magPrev - magNext) * result[0];
        return result;

    }      

}