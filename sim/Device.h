//
// Created by ltb on 25-4-20.
//

#ifndef DEVICE_H
#define DEVICE_H

#endif //DEVICE_H

#include "../function/commom.h"

struct Device {
    vec x,y,z,lambda;

    mat Exin,Eyin;

    vector<mat> Index;
};