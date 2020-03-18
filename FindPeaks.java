
class FindPeaks {

    public static void main(String[] args) {

    }

    public static int[] findPeaks(double[] array) {

        // Only a maximum of N/2 peaks possible (+1 to account for truncation)
        int nProv = (array.length + 1) / 2;
        int[] vAccu = new int[nProv];

        int nPeaks = 0;

        // Edge case: Array has either 1 or 0 entries. -> return 0
        if (array.length < 2) {
            int[] vPeaks = {0};
            return vPeaks;
        }

        // Edge case: Array has 2 entries: 1>2 / 2>1 / 1==2
        if (array.length == 2) {
            if (array[0] > array[1]) {
                vAccu[nPeaks] = 0;
                nPeaks++;
            } else if (array[0] < array[1]) {
                vAccu[nPeaks] = 1;
                nPeaks++;
            } else {
                vAccu[nPeaks] = 0;
                nPeaks++;
                vAccu[nPeaks] = 1;
                nPeaks++;
            }
            return accu2Peaks(vAccu, nPeaks);
        }

        // Avoid cyclic confusion #1
        if (array[0] > array[array.length-1] && array[0] > array[1]) {
            vAccu[nPeaks] = 0;
            nPeaks++;
        }

        // All other cases
        for (int iIdx = 1; iIdx < array.length-1; iIdx++) {

            if (iIdx == 0 && array[iIdx] > array[iIdx+1]) {
                vAccu[nPeaks] = iIdx;
                nPeaks++;
            } else if (array[iIdx] > array[iIdx-1] && array[iIdx] > array[iIdx+1]) {
                vAccu[nPeaks] = iIdx;
                nPeaks++;
            }
        
        }

        if (array[array.length-1] > array[array.length-2] && array[array.length-1] > array[0]) {
            vAccu[nPeaks] = array.length-1;
            nPeaks++;
        }
        
        return accu2Peaks(vAccu, nPeaks);
    }

    public static int[] accu2Peaks(int[] accu, int nPeaks) {

        int[] vPeaks = new int[nPeaks];

        for (int iPeak = 0; iPeak < nPeaks; iPeak++) {
            vPeaks[iPeak] = accu[iPeak];
        }
        return vPeaks;
    }

}