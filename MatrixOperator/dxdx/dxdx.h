//
// Created by ltb on 25-4-20.
//

#ifndef DXDX_H
#define DXDX_H

#endif //DXDX_H

#include "dxdxbase.h"


inline  sp_mat dxdx(int nx , int ny, double dx  ,double dy ) {
    vec I(nx*ny,fill::ones);
    return dxdxbase(I,nx,ny,dx,dy);
}

inline sp_mat dxdx(vec & q , vec & r , int nx , int ny , double dx , double dy) {
    sp_mat R(nx*ny,nx*ny);
    R.diag(0) = r;
    return dxdxbase(q,nx,ny,dx,dy)*R;
}

inline sp_mat dxdx(vec & p , vec &q , vec & r , int nx , int ny , double dx , double dy) {
    int ng = nx * ny ;
    sp_mat P(ng,ng),R(ng,ng);
    P.diag(0)=p;
    R.diag(0)=R;
    return P*dxdxbase(q,nx,ny,dx,dy)*R;
}