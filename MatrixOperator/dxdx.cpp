//
// Created by ltb on 25-4-19.
//
#pragma once
#include "dxdx.h"

#include <ranges>


sp_mat dxdx(double dx, int nx, int ny) {

    const int ng = nx * ny;

    sp_mat S(ng,ng);

    vec a(ng-1,fill::value(1.0 ));
    vec b(ng,fill::value(-2.0));

    // 迪利克雷边界
    for (int i = nx-1;i<ng-1;i+=nx) {
        a(i)=0;
    }

    S.diag(0)=b;
    S.diag(1)=a;
    S.diag(-1)=a;

    return S/dx/dx;
}

sp_mat dxdx(vec q, vec r, double dx, int nx, int ny) {

    const int ng = nx * ny;
    sp_mat S(ng,ng);

    vec z(1);

    vec qr = join_cols(z , q.rows(0,ng-2));
    vec ql = join_cols(q.rows(1,ng-1),z);

    vec qa= 1.0 / ( q.rows(0,ng-2) + q.rows(1,ng-1) );
    for (int i = nx-1;i<ng-1;i+=nx) {
        qa(i)=0;
    }

    vec qb = 1/(ql+q)+1/(q+qr);

    S.diag(0) =- qb;
    S.diag(1)=qa;
    S.diag(-1)=qa;

    sp_mat S1(ng,ng);

    S1.diag(0)=r;

    S = S*S1*2/dx/dx;

    return S;

}
