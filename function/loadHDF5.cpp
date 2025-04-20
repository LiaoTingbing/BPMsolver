//
// Created by ltb on 25-4-20.
//
#pragma once
#include "loadHDF5.h"

void loadData(map<string, cube> &dev,
              string filepath,
              field<string> &s) {
    // map<string,cube> device;
    // device["z"]=aload()
    cube tmp;

    for (int i = 0; i < s.size(); i++) {
        cout << s(i) << endl;
        tmp.load(hdf5_name(filepath, s(i)));
        cout << tmp.n_rows << " " << tmp.n_cols << " " << tmp.n_slices << endl;
        dev[s(i)] = tmp;
        // tmp.print();
    }
}
