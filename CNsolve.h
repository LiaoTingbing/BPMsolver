#pragma once

#include "common.h"


//���ȫʸ����ֹ�ʽ
void  CNFullVectorSolve(DiagStruct& Ax, DiagStruct& Ay,
	DiagStruct& Bx, DiagStruct& By,
	DiagStruct& C, DiagStruct& D,
	cx_vec& u2, cx_vec& v2,
	cx_vec& u1, cx_vec& v1,
    double a, double b);


// ����ʸ����ֹ�ʽ
void CNSemiVectorSolve(cx_vec &u2, 
	cx_vec& u1 ,
	sp_cx_mat& Ax1, sp_cx_mat& Ay1,
	sp_cx_mat& Ax2, sp_cx_mat& Ay2,
	double a, double b);



// ����ϡ������������ĳ˻� A*u
cx_vec sparseMatrixMultipliedByVector(DiagStruct & A ,cx_vec & u);


// ����ϡ������������ĳ˻� ��1+aA)*u
cx_vec sparseMatrixMultipliedByVector(cx_double a , const DiagStruct& A, const cx_vec& u);

// ����	1+a*A     ����ϡ�����
DiagStruct coefficientSparseMatrix(cx_double a , const DiagStruct& A);