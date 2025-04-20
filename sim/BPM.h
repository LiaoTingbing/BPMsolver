//
// Created by ltb on 25-4-20.
//

#ifndef BPM_H
#define BPM_H
#include <map>
#include <string>

#include "../function/commom.h"

class BPM {
    std::map<std::string, cube> *dev;

public:
    BPM(map<string, cube> *dev_);

    BPM();

    ~BPM();

    void init();
};


#endif //BPM_H
