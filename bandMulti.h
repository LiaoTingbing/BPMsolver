#pragma once


#include "common.h"

//计算 y = A*x   A是三对角矩阵，位置为q
// a=-q b=0 c=q
cx_vec triMulti(
	const cx_vec& a, 
	const cx_vec& b, 
	const cx_vec& c
	,int q, 
	const cx_vec& x);

//计算 y = A*x   A是五对角矩阵，位置为q
// a = -q-1 b=-q+1 c=0 d=q-1 e=q+1
cx_vec fiveMulti(const cx_vec & a,
	const cx_vec& b,
	const cx_vec& c,
	const cx_vec& d, 
	const cx_vec& e,
	int q ,
	const cx_vec& x);


