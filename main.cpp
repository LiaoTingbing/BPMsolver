
#include "function/commom.h"
#include "function/loadHDF5.h"
// #include "MatrixOperator/test.h"
#include "sim/BPM.h"

int main() {
    string filePath = "../lumerical/testfile.h5";

    field<string> s = {
        "x", "y", "z", "lambda",
        "index_x", "index_y", "index_z",
        "Exin", "Eyin","neff"
    };

    map<string, cube> dev;
    loadData(dev, filePath, s);

    BPM bpm(&dev);
    bpm.init();
    bpm.propagate();



    // test();



    return 0;
}
