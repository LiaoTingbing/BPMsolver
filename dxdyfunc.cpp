#include "dxdyfunc.h"

sp_cx_mat dxdyfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy)
{

	int ng = nx * ny;

	cx_vec a = p % join_cols(cx_vec(nx + 1), r.rows(0, ng - nx - 2)) % join_cols(cx_vec(1), 1 / q.rows(0, ng - 2));
	cx_vec b = -p % join_cols(cx_vec(nx - 1), r.rows(0, ng - nx)) % join_cols(1.0 / q.rows(1, ng - 1), cx_vec(1));
	cx_vec c = -p % join_cols(r.rows(nx - 1, ng - 1), cx_vec(nx - 1)) % join_cols(cx_vec(1), 1 / q.rows(0, ng - 2));
	cx_vec d = p % join_cols(r.rows(nx + 1, ng - 1), cx_vec(nx + 1)) % join_cols(1 / q.rows(1, ng - 1), cx_vec(1));

	for (int i = 0;i < ng;i += nx) {
		a(i) = 0;
		c(i) = 0;

	}
	for (int i = nx - 1;i < ng;i += nx) {
		b(i) = 0;
		d(i) = 0;
	}


	sp_cx_mat s = spdiags(join_rows(d, c, b, a) / 4 / dx / dy,
		ivec{ -nx - 1, -nx + 1, nx - 1, nx + 1 }, ng, ng).st();

	return s;
}

sp_cx_mat dxdyfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy, DiagStruct& diagv)
{
	int ng = nx * ny;

	cx_vec a = p % join_cols(cx_vec(nx + 1), r.rows(0, ng - nx - 2)) % join_cols(cx_vec(1), 1 / q.rows(0, ng - 2));
	cx_vec b = -p % join_cols(cx_vec(nx - 1), r.rows(0, ng - nx)) % join_cols(1.0 / q.rows(1, ng - 1), cx_vec(1));
	cx_vec c = -p % join_cols(r.rows(nx - 1, ng - 1), cx_vec(nx - 1)) % join_cols(cx_vec(1), 1 / q.rows(0, ng - 2));
	cx_vec d = p % join_cols(r.rows(nx + 1, ng - 1), cx_vec(nx + 1)) % join_cols(1 / q.rows(1, ng - 1), cx_vec(1));

	//for (int i = 0;i < ng;i += nx) {
	//    a(i) = 0;
	//    c(i) = 0;
	//}
	for (int i = nx - 1;i < ng;i += nx) {
		b(i) = 0;
		d(i) = 0;
		a(i - nx + 1) = 0;
		c(i - nx + 1) = 0;
	}


	diagv.diagArr.set_size(ng, 5);
	diagv.diagArr.col(0) = d / 4 / dx / dy;
	diagv.diagArr.col(1) = c / 4 / dx / dy;
	diagv.diagArr.col(2) = cx_vec(ng);
	diagv.diagArr.col(3) = b / 4 / dx / dy;
	diagv.diagArr.col(4) = a / 4 / dx / dy;

	diagv.diagIndex = { -nx - 1, -nx + 1, 0,nx - 1, nx + 1 };

	//sp_cx_mat s = spdiags(diagv.D, diagv.pos, ng, ng).st();

	return sp_cx_mat{};

}

