#pragma once



#include <iostream>
#define ARMA_USE_SUPERLU
#define ARMA_USE_HDF5
#include <armadillo>
#include <ctime>
#include <omp.h>

using namespace std;
using namespace arma;


using namespace  std;
using namespace arma;

const double pi = 3.14159265358979323846;
const double eps0 = 8.854187817e-12;
const double mu0 = 4e-7 * pi;
const double c0 = 1 / sqrt(eps0 * mu0);
const std::complex<double> iu(0, 1);


struct  DiagStruct
{
	// d c b a ÄæÐò Éú³ÉÏ¡Êè¾ØÕó·½±ã
	cx_mat diagArr;
	ivec diagIndex;
};