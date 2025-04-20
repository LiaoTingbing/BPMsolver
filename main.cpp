
// #pragma once

#include "function/commom.h"
#include "function/loadHDF5.h"
#include "MatrixOperator/test.h"
// #include "sim/Device.h"

int main() {

    string filePath = "../lumerical/testfile.h5";

    field<string> s={ "x","y","z","lambda",
        "index_x","index_y","index_z",
    "Exin","Eyin"};

    map<string,cube> dev= loadData( filePath ,  s) ;



    return 0;
}