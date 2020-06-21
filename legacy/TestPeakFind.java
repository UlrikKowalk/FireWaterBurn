
class TestPeakFind {

    public static void main(String[] args) {

        double[] a = {0.2360, -0.8352, -1.2760, 0.6170, 0.6127, 0.2894, 0.3953, -0.8706, -0.4977, -0.1067, -0.6878, 0.3319, 2.3652, -0.4822, 0.6474, -1.0344};

        int[] vPeaks = FindPeaks.findPeaks(a);

        System.out.println("Number of peaks: " + vPeaks.length);

        for (int iPeak = 0; iPeak < vPeaks.length; iPeak++) {
            System.out.println("Peak[" + iPeak + "]: " + vPeaks[iPeak]);
        }
    }

}