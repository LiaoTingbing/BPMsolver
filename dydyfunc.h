#pragma once

#include "common.h"


sp_cx_mat dydyfunc(const cx_vec& p, const  cx_vec& q, const cx_vec& r,
	int nx, int ny, double dx = 1, double dy = 1);

sp_cx_mat dydyfunc(const cx_vec& p, const  cx_vec& q, const cx_vec& r,
	int nx, int ny, double dx , double dy  ,
	DiagStruct & diagv);