//
// Created by ltb on 25-4-20.
//

#ifndef TEST_H
#define TEST_H

#endif //TEST_H


void testdxdybase() {
    // double dx=0.004;
    // double dy=0.004;

    int nx=3;
    int ny = 3;
    int ng = nx * ny ;

    vec q = linspace(1, ng ,ng);

    sp_mat s =  dxdybase(q,nx,ny)  ;
    mat ss = conv_to<mat>::from(s);
    ss.print();
}