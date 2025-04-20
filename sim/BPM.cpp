//
// Created by ltb on 25-4-20.
//

#include "BPM.h"

BPM::BPM(map<string, cube> *dev_) {
    dev = dev_;
}

BPM::BPM() {
}

BPM::~BPM() {
}

void BPM::init() {
    Exin = vectorise((*dev)["Exin"]);
    Eyin = vectorise((*dev)["Eyin"]);
    lambda = (*dev)["lambda"](0);
    neff = (*dev)["neff"](0);
    nx = (*dev)["x"].size();
    ny = (*dev)["y"].size();
    nz = (*dev)["z"].size();
    // Exin.print();
    x = vectorise((*dev)["x"]);
    y = vectorise((*dev)["y"]);
    z = vectorise((*dev)["z"]);
    dx = x(1)-x(0);
    dy = y(1)-y(0);
    dz = z(1)-z(0);
    k0 = 2* pi / lambda;
}

void BPM::propagate() {

}
