#pragma once
#include "Thomas.h"

// Thomas算法函数
vec thomas(const vec& a, const vec& b, const vec& c, const vec& r) {
//     a：下对角线元素（长度为n-1）。
// b：主对角线元素（长度为n）。
// c：上对角线元素（长度为n-1）。
// d：方程组右侧的常数项（长度为n）。
    const size_t n = r.size(); // 方程组的大小
    vec u(n);
    vec g(n);
    double beta = b(0);
    u(0) = r(0)/beta;

    for (int j = 1; j < n; j++) {
        g(j) = c(j-1)/beta;
        beta = b(j) - a(j-1)*g(j);
        u(j) =( r(j) - a(j-1)*u(j-1) ) /beta;
    }

    for (int k = n-2; k>=0; k--) {
        u(k) = u(k) - g(k+1)*u(k+1);
    }
    return u;
}

