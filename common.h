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
const double c0 = 1/sqrt(eps0 * mu0);
const std::complex<double> iu(0, 1);


struct  DiagStruct
{
	// a b c d
	cx_mat diagArr;
	int pos;
};


// 广角系数
const field<cx_vec> Mn{
	{1.0 / 2.0,	0.0 }, //1,0
	{1.0 / 2.0,	0.0 },	//
	{1.0 / 4.0,	1.0 / 2.0,	0.0},
	{3.0 / 32.0,	1.0 / 2.0,	1.0 / 2.0,	0.0},
	{8 / 256.0,		80 / 256.0,		8*24 / 256.0,	8*16 / 256.0,	0 / 256.0},
	{10 / 1024.0,	160 / 1024.0,	672 / 1024.0,	1024 / 1024.0,	512 / 1024.0,	0 / 1024.0}
};
const  field<cx_vec> Nn{
	{0.0,			1.0},
	{1.0/4.0,		1.0},
	{1.0/16.0,	3.0/4.0,	1.0},
	{1.0/64.0,	3.0/8.0,	5.0/4.0,		1.0},
	{1 / 256.0,		40 / 256.0,		240 / 256.0,	448 / 256.0,	256 / 256.0},
	{1 / 1024.0,	60 / 1024.0,	560 / 1024.0,	1792 / 1024.0,	2304 / 1024.0,	1024/ 1024.0}
};

