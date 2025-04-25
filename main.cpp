// bpm20250421.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "common.h"
#include "loadHDF5.h"
#include "BPM.h"
#include "dxdxfunc.h"
#include "dydyfunc.h"

//void test() {
//    vec I(16, fill::ones);
//    sp_mat s = dydyfunc(I, I, I, 4,4);
//    mat ss = conv_to<mat>::from(s);
//    ss.print();
//
//}

int main()
{
    //omp_set_num_threads(6);
    //test();

    string filePathrosft = "lumerical/rsoft.h5";


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
    bpm.FullVector_propagate();
    bpm.postData();
    //bpm.getMatrix();
 
    //bpm.Qusi_TM_Propagate();    //  EX
    //bpm.Qusi_TE_Propagate();    //  EY
    //bpm.postData();
    //bpm.compute_FullVectorMatrix();

}

