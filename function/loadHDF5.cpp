//
// Created by ltb on 25-4-20.
//

#include "loadHDF5.h"

map<string,cube> loadData(string filepath, field<string> &s) {

    map<string,cube> device;
    // device["z"]=aload()
    cube tmp ;

    for ( int i = 0 ;i < s.size();i++) {
        cout<<s(i) <<endl;
        tmp.load(hdf5_name(filepath,s(i)));
        cout<<tmp.n_rows<<" " <<tmp.n_cols<<" "<<tmp.n_slices<<endl;
        device[s(i)] = tmp;
        // tmp.print();
    }

    return  device;
}
