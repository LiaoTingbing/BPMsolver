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
    Exin = vectorise((*dev)["Exin"])  + 0.0* iu;
    Eyin = vectorise((*dev)["Eyin"]) + 0.0* iu;
    Et = join_cols(Exin, Eyin) + 0.0* iu;
    lambda = (*dev)["lambda"](0);
    n0 = (*dev)["neff"](0);
    nx = (*dev)["x"].size();
    ny = (*dev)["y"].size();
    nz = (*dev)["z"].size();
    // Exin.print();
    x = vectorise((*dev)["x"]);
    y = vectorise((*dev)["y"]);
    z = vectorise((*dev)["z"]);
    dx = x(1) - x(0);
    dy = y(1) - y(0);
    dz = z(1) - z(0);
    k0 = 2 * pi / lambda;
}




void BPM::propagate() {


    int i =0 ;

    vec ERX =pow( vectorise (  (*dev)["index_x"].slice(i)  )  ,2);
    vec ERY =pow( vectorise (  (*dev)["index_y"].slice(i)  ) ,2);
    vec ERZ =pow(  vectorise (  (*dev)["index_z"].slice(i)  ) ,2);
    vec ERX2 =pow( vectorise (  (*dev)["index_x"].slice(i+1)  )  ,2);
    vec ERY2 =pow( vectorise (  (*dev)["index_y"].slice(i+1)  ) ,2);
    vec ERZ2 =pow(  vectorise (  (*dev)["index_z"].slice(i+1)  ) ,2);

    vec ziv(nx*ny);

    onestepViolence(ERX, ERY, ERZ, ziv, ziv, ERX2, ERY2, ERZ2,ziv,
                          ziv,  Exin,  Eyin,   Et,   dx,   dy,   dz,   k0,
                           n0,
                           alpha,   nx,   ny,
                           nz);



}
