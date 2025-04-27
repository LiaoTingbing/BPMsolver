#pragma once

#include "common.h"
#include "cn_solve.h"
 

#include "dxdxfunc.h"
#include "dydyfunc.h"
#include "dydxfunc.h"
#include "dxdyfunc.h"

#include "thomas.h"



class Bpm
{
	map<string, cube>* dev_;

	double dx_, dy_, dz_;
	int nx_, ny_, nz_,nt_;
	double lambda_, k0_, n0_;
	double alpha_;              // CN 差分控制参数
 
	field<DiagStruct> Ax_;           //依据书上命名
	field<DiagStruct> Ay_;
	field<DiagStruct> Bx_;
	field<DiagStruct> By_;
	field<DiagStruct> C_;
	field<DiagStruct> D_;


	cx_vec sx_,sy_,isx_,isy_;
 

	//	首字母大写矩阵
	//	首字母小写列向量
	cx_mat ex_;      //  nt*nz
	cx_mat ey_;		//	nt *nz

public:

	Bpm(map<string, cube>* dev) {
		dev_ = dev;
	};

	Bpm() {
	};

	~Bpm() {
	};

	
	void init();

	void computePML(int layersPML = 10);      // 计算PML参数

	void computeMatrix();
 
	void qusiTmPropagate();

	void qusiTePropagate();

	void fullVectorPropagateSimple();

	// Order 0 = 1,0;
	// Order 1 = 1,1;
	// Order 2 = 2,2;
	// Order 3 = 3,3;
	// Order 4 = 4,4;
	// Order 5 = 5,5;
	void fullVectorWideAnglePropagateSimple(int order=0);

	void postData();

};

