#include "dxdxfunc.h"

sp_cx_mat dxdxfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r,
	int nx, int ny,
	double dx, double dy)
{
	int ng = nx * ny;
	cx_vec a = p % join_cols(cx_vec(1), r.rows(0, ng - 2))
		/ (join_cols(cx_vec(1), q.rows(0, ng - 2)) + q);

	cx_vec b = -p % r % (
		1 / (join_cols(cx_vec(1), q.rows(0, ng - 2)) + q)
		+ 1 / (q + join_cols(q.rows(1, ng - 1), cx_vec(1)))
		);
	cx_vec c = p % join_cols(r.rows(1, ng - 1), cx_vec(1))
		/ (q + join_cols(q.rows(1, ng - 1), cx_vec(1)));

	for (int i = 0;i < ng;i += nx)
		a(i) = 0;

	for (int j = nx - 1;j < ng;j += nx)
		c(j) = 0;

	sp_cx_mat s = spdiags(join_rows(c, b, a) * 2 / dx / dx,
		ivec{ -1,0,1 }, ng, ng).st();
	return s;

}

sp_cx_mat dxdxfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy, DiagStruct& diagv)
{
	int ng = nx * ny;
	cx_vec a = p % join_cols(cx_vec(1), r.rows(0, ng - 2))
		/ (join_cols(cx_vec(1), q.rows(0, ng - 2)) + q);

	cx_vec b = -p % r % (
		1 / (join_cols(cx_vec(1), q.rows(0, ng - 2)) + q)
		+ 1 / (q + join_cols(q.rows(1, ng - 1), cx_vec(1)))
		);
	cx_vec c = p % join_cols(r.rows(1, ng - 1), cx_vec(1))
		/ (q + join_cols(q.rows(1, ng - 1), cx_vec(1)));

	//for (int i = 0;i < ng;i += nx)
	//	a(i) = 0;


	//for (int j = nx - 1;j < ng;j += nx)
	//	c(j) = 0;

	for (int j = nx - 1;j < ng;j += nx)
	{
		c(j) = 0;
		a(j - nx + 1) = 0;
	}

	diagv.diagArr.set_size(ng, 3);
	diagv.diagArr.col(0) = c * 2 / dx / dx;
	diagv.diagArr.col(1) = b * 2 / dx / dx;
	diagv.diagArr.col(2) = a * 2 / dx / dx;

	diagv.diagIndex = { -1 ,0,1 };

	//sp_cx_mat s = spdiags(diagv.D,diagv.pos, ng, ng).st();
	return sp_cx_mat{};
}

