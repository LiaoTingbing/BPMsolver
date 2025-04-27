#pragma once

#include "common.h"

// Thomas sove 三对角
// a=-1 b=0 c=1
cx_vec thomasAlgorithm(
    const cx_vec& a, const cx_vec& b, const cx_vec& c,
    const cx_vec& d);


//  托马斯求解   a=-p , b=0, c=p
cx_vec thomasAlgorithm(
    const cx_vec& a, const cx_vec& b, const cx_vec& c,
    int p, const cx_vec& d);


 