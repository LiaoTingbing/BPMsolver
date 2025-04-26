#include "CNsolve.h"

cx_vec thomasSolve(const DiagStruct& A, const cx_vec& y)
{
    
    return thomas_algorithm(
        A.diagArr.col(2),   //a
        A.diagArr.col(1),   //b
        A.diagArr.col(0),   //c
        A.pos,              // 对角位置
        y);
}

cx_vec sparseMatrixMultipliedByVector(const DiagStruct& D, const cx_vec& u)
{
    // 计算稀疏矩阵与向量的乘积 D*u
    return fiveMulti(
        D.diagArr.col(4),//a
        D.diagArr.col(3),//b
        D.diagArr.col(2),//c
        D.diagArr.col(1),//d
        D.diagArr.col(0),//e
        D.pos,
        u);
}

cx_vec sparseMatrixMultipliedByVector(cx_double a, const DiagStruct& A, const cx_vec& u)
{
    // 计算稀疏矩阵与向量的乘积 （1+aA)*u
    return triMulti(
		a * A.diagArr.col(2),//a
		1 + a * A.diagArr.col(1),//b
		a * A.diagArr.col(0),//c
		A.pos,
		u);
}

DiagStruct coefficientSparseMatrix(cx_double a, const DiagStruct& A)
{
    return   {
        join_rows(a * A.diagArr.col(0),1 + a * A.diagArr.col(1),a * A.diagArr.col(2)),
        A.diagIndex,
        A.pos
    };
}

cx_vec CNfisrtOne(const DiagStruct& Ayr, const DiagStruct& Ayl, const cx_vec& ur, cx_double a, cx_double b)
{
    // 计算(1+aAyr)*ur
    return thomasSolve(
        coefficientSparseMatrix(b, Ayl),
        sparseMatrixMultipliedByVector(a, Ayr, ur)
    )  ;
}

