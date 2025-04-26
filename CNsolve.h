#pragma once

#include "common.h"
#include "bandMulti.h"
#include "thomas.h"


//	求解三对角矩阵方程 Ax = b Thomas
// 输入 A是三对角矩阵
cx_vec thomasSolve( const DiagStruct & A , const cx_vec & y);

// 计算稀疏矩阵与向量的乘积 D*u
// D是五对角
cx_vec sparseMatrixMultipliedByVector(const DiagStruct& D, cx_vec& u);


// 计算稀疏矩阵与向量的乘积 （1+aA)*u
cx_vec sparseMatrixMultipliedByVector(cx_double a, const DiagStruct& A, const cx_vec& u);

// 计算	1+a*A     返回稀疏矩阵
DiagStruct coefficientSparseMatrix(cx_double a, const DiagStruct& A);


// CN差分第一步1.1
// 计算（1+bAyl)*ul = (1+aAyr)*ur
cx_vec CNfisrtOne(const DiagStruct& Ayr, const DiagStruct& Ayl,
	const cx_vec& ur, cx_double a, cx_double b);

// CN差分第一步1.2
cx_vec CNfisrtTwo(
	cx_double a, 
	const DiagStruct& Byr, const DiagStruct& Dr,
	const  cx_vec& vr, const cx_vec& ur,
	const cx_double b,
	const DiagStruct Byl,const DiagStruct& Dyl,
	const cx_vec ul);