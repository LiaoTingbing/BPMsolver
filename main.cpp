// bpm20250421.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "common.h"
#include "loadHDF5.h"
#include "BPM.h"


int main()
{
    string filePath = "lumerical/testfile.h5";

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
    bpm.postData();

}

