#pragma once

#include "common.h"
#include "bandMulti.h"


//求解全矢量差分公式
void  CNFullVectorSolve(DiagStruct& Ax, DiagStruct& Ay,
	DiagStruct& Bx, DiagStruct& By,
	DiagStruct& C, DiagStruct& D,
	cx_vec& u2, cx_vec& v2,
	cx_vec& u1, cx_vec& v1,
    double a, double b);


// 求解半矢量差分公式
void CNSemiVectorSolve(cx_vec &u2, 
	cx_vec& u1 ,
	sp_cx_mat& Ax1, sp_cx_mat& Ay1,
	sp_cx_mat& Ax2, sp_cx_mat& Ay2,
	double a, double b);



// 计算稀疏矩阵与向量的乘积 D*u
// D是五对角
cx_vec sparseMatrixMultipliedByVector(DiagStruct & D ,cx_vec & u);


// 计算稀疏矩阵与向量的乘积 （1+aA)*u
cx_vec sparseMatrixMultipliedByVector(cx_double a , const DiagStruct& A, const cx_vec& u);

// 计算	1+a*A     返回稀疏矩阵
DiagStruct coefficientSparseMatrix(cx_double a , const DiagStruct& A);


