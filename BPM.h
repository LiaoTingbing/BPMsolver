#pragma once

#include "common.h"
#include "CNsolve.h"
 

#include "dxdxfunc.h"
#include "dydyfunc.h"
#include "dydxfunc.h"
#include "dxdyfunc.h"

#include "thomas.h"



class BPM
{
	map<string, cube>* dev;

	double dx, dy, dz;
	int nx, ny, nz,nt;
	double lambda, k0, n0;
	double alpha;              // CN 差分控制参数

	//field<sp_cx_mat> Ax;
	//field<sp_cx_mat> Ay;
	//field<sp_cx_mat> Bx;
	//field<sp_cx_mat> By;
	//field<sp_cx_mat> C;
	//field<sp_cx_mat> D;

	field<DiagStruct> Ax_V;
	field<DiagStruct> Ay_V;
	field<DiagStruct> Bx_V;
	field<DiagStruct> By_V;
	field<DiagStruct> C_V;
	field<DiagStruct> D_V;


	cx_vec sx,sy,isx,isy;

	int layersPML;


 
	//vec exin;
	//vec eyin;
 

	//	首字母大写矩阵
	//	首字母小写列向量
	cx_mat Ex;      //  nt*nz
	cx_mat Ey;		//	nt *nz

public:

	BPM(map<string, cube>* dev_) {
		dev = dev_;
	};
	BPM() {

	};
	~BPM() {

	};

	
	void init();


	void compute_PML();      // 计算PML参数

	void compute_Matrix();
 
	void Qusi_TM_Propagate();

	void Qusi_TE_Propagate();

	void FullVector_propagate();

	void FullVector_propagate_simple();


	void postData();


 
	
};

