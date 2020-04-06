package com.JadeHS.Viwer;

import java.util.*;

public class CandidateList {

    private ArrayList<ArrayList<Double>> list;

    public CandidateList() {
        this.list = new ArrayList<>();
    }

    public void add(ArrayList<Double> entry) {
        this.list.add(entry);
    }

    public ArrayList<Double> get(int idx) {
        return this.list.get(idx);
    }

    public int size() {
        return this.list.size();
    }
}