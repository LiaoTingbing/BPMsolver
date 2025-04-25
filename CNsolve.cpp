#include "CNsolve.h"

void CNFullVectorSolve(cx_vec& u2, cx_vec& v2, cx_vec& u1, cx_vec& v1, sp_cx_mat& Ax1, sp_cx_mat& Ay1, sp_cx_mat& Bx1, sp_cx_mat& By1, sp_cx_mat& C1, sp_cx_mat& D1, sp_cx_mat& Ax2, sp_cx_mat& Ay2, sp_cx_mat& Bx2, sp_cx_mat& By2, sp_cx_mat& C2, sp_cx_mat& D2, double a, double b)
{
   sp_mat II = speye(size(Ax1));

    //% 第一步
       cx_vec utmp = spsolve(II + b * Ay2, 
            (II + a * Ay1) * u1);

    cx_vec vtmp = spsolve(II + b * By2,  
        (II + a * By1) * v1 + a * D1 * u1 - b * D2 * utmp);

    //% 第二步
        v2 = spsolve(II + b * Bx2,  
            (II + a * Bx1) * vtmp);
    u2 = spsolve(II + b * Ax2,  
        (II + a * Ax1) * utmp + a * C1 * vtmp - b * C2 * v2);
}

void CNSemiVectorSolve(cx_vec& u2, cx_vec& u1, sp_cx_mat& Ax1, sp_cx_mat& Ay1, sp_cx_mat& Ax2, sp_cx_mat& Ay2, double a, double b)
{
    sp_mat II = speye(size(Ax1));
    cx_vec utmp = spsolve(II + b * Ay2, (II + b * Ay2)*u1);
	u2 = spsolve(II + b * Ax2, (II + a * Ax1) * utmp);
}

