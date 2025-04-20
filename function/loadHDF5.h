//
// Created by ltb on 25-4-20.
//

#ifndef LOADHDF5_H
#define LOADHDF5_H

#pragma once
#include "commom.h"



void loadData(map<string, cube> &dev,
              string filepath,
              field<string> &s);



#endif //LOADHDF5_H



