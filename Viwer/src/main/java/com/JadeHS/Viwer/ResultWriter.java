package com.JadeHS.Viwer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class ResultWriter {

    private String filename;
    private File file;

    public ResultWriter(String filename) {
        this.filename = filename;
        this.file = new File(this.filename);
        this.delete();
    }

    public void delete() {
        this.file.delete();
    }

    public boolean write(double data) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(this.filename, true));
            writer.append(data + " ");
            writer.close();
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    public boolean write(double[] data) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(this.filename, true));
            for (int ival = 0; ival < data.length; ival++) {
                writer.append(data[ival] + " ");
            }
            writer.close();
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    public boolean write(String data) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(this.filename, true));
            writer.append(data);
            writer.close();
            return true;
        } catch (Exception e) {
            return false;
        }
    }

}