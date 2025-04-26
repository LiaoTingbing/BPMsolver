#include "CNsolve.h"

 

void CNSemiVectorSolve(cx_vec& u2, cx_vec& u1, sp_cx_mat& Ax1, sp_cx_mat& Ay1, sp_cx_mat& Ax2, sp_cx_mat& Ay2, double a, double b)
{
    sp_mat II = speye(size(Ax1));
    cx_vec utmp = spsolve(II + b * Ay2, (II + b * Ay2)*u1);
	u2 = spsolve(II + b * Ax2, (II + a * Ax1) * utmp);
}

cx_vec sparseMatrixMultipliedByVector(DiagStruct& A, cx_vec& u)
{
    // 计算稀疏矩阵与向量的乘积 A*u
    return spdiags(A.diagArr, A.diagIndex, 
        A.diagArr.n_rows, A.diagArr.n_rows).st()* u;
}

cx_vec sparseMatrixMultipliedByVector(cx_double a, const DiagStruct& A, const cx_vec& u)
{
    // 计算稀疏矩阵与向量的乘积 （1+aA)*u
    return spdiags( 
        join_rows(a*A.diagArr.col(0),1+ a * A.diagArr.col(1), a * A.diagArr.col(2) )
        , A.diagIndex,
        A.diagArr.n_rows, A.diagArr.n_rows).st() * u;
}

DiagStruct coefficientSparseMatrix(cx_double a, const DiagStruct& A)
{
    return   {
        join_rows(a * A.diagArr.col(0),1 + a * A.diagArr.col(1),a * A.diagArr.col(2)),
        A.diagIndex
    };
}

