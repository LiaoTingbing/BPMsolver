#pragma once

#include "common.h"


sp_cx_mat dxdxfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r,
	int nx, int ny,
	double dx, double dy );

sp_cx_mat dxdxfunc(const cx_vec& p, const cx_vec& q, const cx_vec& r,
	int nx, int ny,
	double dx, double dy  , DiagVector & diagv);
 