

class TestClass {

    public static void main(String[] args) {

        String filename = "media/Boomerang.wav";
        OpenWaveFile wave = new OpenWaveFile(filename);
        wave.open();

        print("Channels: " + wave.getNumChannels() + ", Length: " + wave.getLength() + ", Samplerate: " + wave.getSampleRate());
        
        int blocksize = 512;
        double overlap = 0.0;
        
        long startTime = System.nanoTime();


        Viwer viwer = new Viwer(wave.getSampleRate(), blocksize, overlap, wave.getNumChannels());
        viwer.analyseAndFilter(wave.returnAudio());
   

        long endTime = System.nanoTime();
        long duration = (endTime - startTime); 
        print("Execution took " + duration/1e6 + "ms");
    }

    public static void print(String string) {
        System.out.println(string);
    }

}