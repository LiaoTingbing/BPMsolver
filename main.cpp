// bpm20250421.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "common.h"
#include "load_hdf5.h"
#include "bpm.h"



int main()
{
    //omp_set_num_threads(6);
    //test();

    string filePathRosft = "matlab/rsoft.h5";


    field<string> rosft = {
    "x", "y", "z", 
	"Epsx", "Epsy", "Epsz","Epsxy","Epsyx",
    "Exin", "Eyin",
    "neff","lambda"
    };

    map<string, cube> dev;
 
    loadData(dev, filePathRosft, rosft);
 

	Bpm bpm(&dev);

    bpm.init();
    bpm.computePML();
    bpm.computeMatrix();
    //bpm.FullVector_propagate_simple();
    //Pade 1，1最优
    bpm.fullVectorWideAnglePropagateSimple(); // Pade 5,5
    bpm.postData();
 

}

