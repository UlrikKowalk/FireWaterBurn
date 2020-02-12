

class TestClass {

    public static void main(String[] args) {

        String filename = "D:/Java/media/45_8_Milbe.wav";
        OpenWaveFile wave = new OpenWaveFile(filename);
        wave.open();

        print("Channels: " + wave.getNumChannels() + ", Length: " + wave.getLength() + ", Samplerate: " + wave.getSampleRate());
        
        int blocksize = 1024;
        double overlap = 0.0;


        Localisation localisation = new Localisation();

        int nRuns = 1;

        long startTime = System.nanoTime();

        for (int iRun = 0; iRun < nRuns; iRun++) {
            localisation.performExperiment(wave.returnAudio(), wave.getSampleRate(), blocksize, overlap);
        }

        long endTime = System.nanoTime();

        long duration = (endTime - startTime)/nRuns; 
        print("Execution took " + duration/1e6 + "ms");




    }

    public static void print(String string) {
        System.out.println(string);
    }

}