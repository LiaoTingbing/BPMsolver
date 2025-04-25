#include "dydyfunc.h"

sp_cx_mat dydyfunc(const cx_vec& p, const  cx_vec& q, const cx_vec& r,
	int nx, int ny, double dx, double dy)
{
	int ng = nx * ny;
 
	cx_vec a = p % join_cols(cx_vec(nx),r.rows(0,ng - nx-1))
		/ ( join_cols(cx_vec(nx),q.rows(0,ng - nx-1)) + q);
	cx_vec b = -p % r % (1 / (join_cols(cx_vec(nx),q.rows(0,ng - nx-1)) + q)
		+ 1 / (q + join_cols(q.rows(nx ,ng-1),cx_vec(nx))));
	cx_vec c = p % join_cols(r.rows(nx,ng-1),cx_vec(nx))
		/ (q + join_cols(q.rows(nx,ng-1),cx_vec(nx)));

	sp_cx_mat s =  spdiags(join_rows(c, b, a) * 2 / dy / dy, 
		ivec{ -nx,0,nx }, ng, ng).st();
	return s;

	  
 
}

sp_cx_mat dydyfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy, DiagVector& diagv)
{
	int ng = nx * ny;

	cx_vec a = p % join_cols(cx_vec(nx), r.rows(0, ng - nx - 1))
		/ (join_cols(cx_vec(nx), q.rows(0, ng - nx - 1)) + q);
	cx_vec b = -p % r % (1 / (join_cols(cx_vec(nx), q.rows(0, ng - nx - 1)) + q)
		+ 1 / (q + join_cols(q.rows(nx, ng - 1), cx_vec(nx))));
	cx_vec c = p % join_cols(r.rows(nx, ng - 1), cx_vec(nx))
		/ (q + join_cols(q.rows(nx, ng - 1), cx_vec(nx)));

 
	diagv.diagArr.set_size(ng, 3);
	diagv.diagArr.col(0) = c * 2 / dy / dy;
	diagv.diagArr.col(1) = b * 2 / dy / dy;
	diagv.diagArr.col(2) = a * 2 / dy / dy;

	diagv.diagIndex = { -nx,0,nx };

	//sp_cx_mat s = spdiags(diagv.D,diagv.pos, ng, ng).st();
	return sp_cx_mat{};
}

 