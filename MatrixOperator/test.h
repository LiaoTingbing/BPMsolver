//
// Created by ltb on 25-4-20.
//

#ifndef TEST_H
#define TEST_H

#endif //TEST_H
#include  "dxdy/dxdybase.h"
#include "dydx/dydxbase.h"
#include "dxdx/dxdxbase.h"
#include "dydy/dydybase.h"
inline  void test() {
    int nx=3;
    int ny = 3;
    int ng = nx * ny ;

    vec q = linspace(1, ng ,ng);

    sp_mat s =  dydxbase(q,nx,ny)  ;
    mat ss = conv_to<mat>::from(s);
    ss.print();

}