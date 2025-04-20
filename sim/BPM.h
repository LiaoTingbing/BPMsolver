//
// Created by ltb on 25-4-20.
//

#ifndef BPM_H
#define BPM_H

#pragma once
#include <map>
#include <string>
#include "../function/commom.h"
#include "../MatrixOperator/matrixHeader.h"

class BPM {
    std::map<std::string, cube> *dev;

    vec x,y,z;
    vec Exin,Eyin;
    double lambda,neff,dx,dy,dz,k0;
    int nx , ny , nz;

public:
    BPM(map<string, cube> *dev_);

    BPM();

    ~BPM();

    void init();
    void propagate();
};


#endif //BPM_H
