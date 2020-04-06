package com.JadeHS.Viwer;

import java.io.File;
import com.JadeHS.Viwer.WavFile.WavFile;


class OpenWaveFile {

    private int samplingrate;
    private String filename;
    private double[][] buffer;
    private int numChannels;
    private int length;
    
    public OpenWaveFile(String filename) {
        this.filename = filename;
    }

    public void open() {
        try
        {
            // Open the wav file specified as the first argument
            WavFile wavFile = WavFile.openWavFile(new File(filename));

            // Get the number of audio channels in the wav file
            this.numChannels = wavFile.getNumChannels();
            this.samplingrate = (int) wavFile.getSampleRate();
            this.length = (int) wavFile.getNumFrames();
            this.buffer = new double[this.numChannels][this.length];

            wavFile.readFrames(buffer, 0, this.length);
            // Close the wavFile
            wavFile.close();

        } catch (Exception e) {e.printStackTrace();}
    }

    public int getSampleRate() {
        return this.samplingrate;
    }

    public double[][] returnAudio() {
        return this.buffer;
    }

    public int getNumChannels() {
        return this.numChannels;
    }

    public int getLength() {
        return this.length;
    }

}