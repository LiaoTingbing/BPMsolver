
// #pragma once

#include "function/commom.h"

#include "MatrixOperator/test.h"

int main() {

    string filePath = "lumerical/testfile.h5";

    vec x;
    x.load(hdf5_name(filePath,"x"));
    x.print();

    // test();


    return 0;
}