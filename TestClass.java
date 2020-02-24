

class TestClass {

    public static void main(String[] args) {

        String filename = "D:/Java/media/0_recording_joined.wav";
        OpenWaveFile wave = new OpenWaveFile(filename);
        wave.open();

        print("Channels: " + wave.getNumChannels() + ", Length: " + wave.getLength() + ", Samplerate: " + wave.getSampleRate());
        
        int blocksize = 512;
        double overlap = 0.0;


        Localisation localisation = new Localisation();

        

        int nRuns = 1;

        long startTime = System.nanoTime();

        //for (int iRun = 0; iRun < nRuns; iRun++) {
            double[][] vResult = localisation.performExperiment(wave.returnAudio(), wave.getSampleRate(), blocksize, overlap);
        //}

        long endTime = System.nanoTime();

        long duration = (endTime - startTime)/nRuns; 
        
        /*for (int iRes = 0; iRes < vResult.length; iRes++) {
            print(" " + vResult[iRes] / Math.PI * 180f);
        }*/
        
        print("Execution took " + duration/1e6 + "ms");


        


    }

    public static void print(String string) {
        System.out.println(string);
    }

}