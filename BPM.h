#pragma once

#include "common.h"
#include "dxdx.h"
#include "dxdy.h"
#include "dydx.h"
#include "dydy.h"


class BPM
{
	map<string, cube>* dev;

	double dx, dy, dz;
	int nx, ny, nz,nt;
	double lambda, k0, n0;
	double alpha;              // CN 差分控制参数

	vec Ein;
	//field<cx_vec> Eout;
	cx_mat Eout;

public:

	BPM(map<string, cube>* dev_) {
		dev = dev_;
	};
	BPM() {

	};
	~BPM() {

	};
	
	void init();

	sp_mat Calculate_P(int i);
 
	void propagate();

	void postData();
	
};

