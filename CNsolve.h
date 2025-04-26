#pragma once

#include "common.h"


//求解全矢量差分公式
void  CNFullVectorSolve(cx_vec& u2, cx_vec& v2, 
    cx_vec & u1, cx_vec& v1,
    sp_cx_mat & Ax1, sp_cx_mat&  Ay1, 
    sp_cx_mat& Bx1, sp_cx_mat& By1, 
    sp_cx_mat& C1, sp_cx_mat& D1,
    sp_cx_mat& Ax2, sp_cx_mat& Ay2,
    sp_cx_mat& Bx2, sp_cx_mat& By2,
    sp_cx_mat& C2, sp_cx_mat& D2,
    double a, double b);


// 求解半矢量差分公式
void CNSemiVectorSolve(cx_vec &u2, 
	cx_vec& u1 ,
	sp_cx_mat& Ax1, sp_cx_mat& Ay1,
	sp_cx_mat& Ax2, sp_cx_mat& Ay2,
	double a, double b);



// 计算稀疏矩阵与向量的乘积 A*u
cx_vec sparseMatrixMultipliedByVector(DiagStruct & A ,cx_vec & u);


// 计算稀疏矩阵与向量的乘积 （1+aA)*u
cx_vec sparseMatrixMultipliedByVector(double a , DiagStruct& A, cx_vec& u);
