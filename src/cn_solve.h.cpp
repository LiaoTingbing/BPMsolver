#include "../include/cn_solve.h"

cx_vec thomasSolve(const DiagStruct& A, const cx_vec& y)
{

	return thomasAlgorithm(
		A.diagArr.col(0),   //a
		A.diagArr.col(1),   //b
		A.diagArr.col(2),   //c
		A.pos,              // 对角位置
		y);
}

cx_vec sparseMatrixMultipliedByVector(const DiagStruct& D, const cx_vec& u)
{
	// 计算稀疏矩阵与向量的乘积 D*u
	return fiveMulti(
		D.diagArr.col(0),//a
		D.diagArr.col(1),//b
		D.diagArr.col(2),//c
		D.diagArr.col(3),//d
		D.diagArr.col(4),//e
		D.pos,
		u);
}

cx_vec sparseMatrixMultipliedByVector(cx_double a, const DiagStruct& A, const cx_vec& u)
{
	// 计算稀疏矩阵与向量的乘积 （1+aA)*u
	return triMulti(
		a * A.diagArr.col(0),//a
		1 + a * A.diagArr.col(1),//b
		a * A.diagArr.col(2),//c
		A.pos,
		u);
}

DiagStruct coefficientSparseMatrix(cx_double a, const DiagStruct& A)
{
	return   {
		join_rows(a * A.diagArr.col(0),1.0 + a * A.diagArr.col(1),a * A.diagArr.col(2)),
		A.pos
	};
}


cx_vec cnFirstOne(cx_double a, const DiagStruct& Ayr, const cx_vec& ur, cx_double b, const DiagStruct Ayl)
{
	return thomasSolve(
		coefficientSparseMatrix(b, Ayl),
		sparseMatrixMultipliedByVector(a, Ayr, ur)
	);
}

cx_vec cnFirstTwo(cx_double a, const DiagStruct& Byr, const DiagStruct& Dr, const cx_vec& vr, const cx_vec& ur, const cx_double b, const DiagStruct Byl, const DiagStruct& Dl, const cx_vec ul)
{
	return thomasSolve(
		coefficientSparseMatrix(b, Byl),
		sparseMatrixMultipliedByVector(a, Byr, vr)
		+ a * sparseMatrixMultipliedByVector(Dr, ur)
		- b * sparseMatrixMultipliedByVector(Dl, ul)
	);
}



cx_vec cnSecondOne(cx_double a, const DiagStruct& Bxr, const cx_vec& vr, const cx_double b, const DiagStruct Bxl)
{
	return thomasSolve(
		coefficientSparseMatrix(b, Bxl),
		sparseMatrixMultipliedByVector(a, Bxr, vr)
	);
}

cx_vec cnSecondTwo(cx_double a, const DiagStruct& Axr, const DiagStruct& Cr, const cx_vec& ur, const cx_vec& vr, const cx_double b, const DiagStruct Axl, const DiagStruct& Cl, const cx_vec vl)
{
	return thomasSolve(
		coefficientSparseMatrix(b, Axl),
		sparseMatrixMultipliedByVector(a, Axr, ur)
		+ a * sparseMatrixMultipliedByVector(Cr, vr)
		- b * sparseMatrixMultipliedByVector(Cl, vl)
	);
}

void cnSolve(cx_double a, cx_double b, const DiagStruct& Ayr, const DiagStruct& Ayl, const DiagStruct& Byr, const DiagStruct& Byl, const DiagStruct& Axr, const DiagStruct& Axl, const DiagStruct& Bxr, const DiagStruct& Bxl, const DiagStruct& Cr, const DiagStruct& Cl, const DiagStruct& Dr, const DiagStruct& Dl, const cx_vec& uin, const cx_vec& vin, cx_vec& uout, cx_vec& vout)
{
	//CN差分第一步1.1
	cx_vec utmp = cnFirstOne(a, Ayr, uin, b, Ayl);
	//CN差分第一步1.2
	cx_vec vtmp = cnFirstTwo(a, Byr, Dr, vin, uin, b, Byl, Dl, utmp);
	//CN差分第二步2.1
	vout = cnSecondOne(a, Bxr, vtmp, b, Bxl);
	//CN差分第二步2.2
	uout = cnSecondTwo(a, Axr, Cr, utmp, vtmp, b, Axl, Cl, vout);
	// 计算结果
}


