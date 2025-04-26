// bpm20250421.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "common.h"
#include "loadHDF5.h"
#include "BPM.h"



int main()
{
    //omp_set_num_threads(6);
    //test();

    string filePathrosft = "matlab/rsoft.h5";


    field<string> rosft = {
    "x", "y", "z", 
	"Epsx", "Epsy", "Epsz","Epsxy","Epsyx",
    "Exin", "Eyin",
    "neff","lambda"
    };

    map<string, cube> dev;
 
    loadData(dev, filePathrosft, rosft);
 

	BPM bpm(&dev);

    bpm.init();
    bpm.compute_PML();
    bpm.compute_Matrix();
    //bpm.FullVector_propagate_simple();
    bpm.FullVector_WideAngle_propagate_simple(1); // Pade 5,5
    bpm.postData();
 

}

