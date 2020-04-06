package com.JadeHS.Viwer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class MatReader {

    private File file;

    public MatReader(File file) {
        this.file = file;
    }

    public void readAll() {
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(this.file));
            String test = bufferedReader.readLine();
            bufferedReader.close();
            System.out.println("result:  " + test);
        } catch (Exception e) {

            
            
        }
        int a = 1;
    }

}