//
// Created by ltb on 25-4-19.
//

#ifndef COMMOM_H
#define COMMOM_H

#endif //COMMOM_H

// #include <hdf5.h>
#define ARMA_USE_HDF5
#pragma once
#include <armadillo>
#include <cstring>
#include <iostream>

using namespace  std;
using namespace arma;

const double pi = 3.14159265358979323846;
const double eps0 = 8.854187817e-12;
const double mu0 = 4e-7*pi;
const double c0 = 1/sqrt(eps0*mu0);
const cx_double iu(0, 1);


