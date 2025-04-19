//
// Created by ltb on 25-4-19.
//

#include "dxdy.h"


// q r point
sp_mat dxdybase(vec &q, int nx, int ny) {
    int ng = nx * ny ;
    sp_mat s(ng ,ng );

    vec z(1);

    vec d= 1 /q(1,0 , size(ng-(nx+1),1));
    for ( int i = nx-1;i<d.size();i+=nx) {
        d(i)=0;
    }

    vec c=-1/ join_cols( z , q(0,0 , size(ng-(nx-1)-1,1))) ;
    for ( int i = 0;i<c.size();i+=nx) {
        c(i)=0;
    }

    vec b =-1/join_cols(z ,  q(1,0,size(ng-(nx-1)-1,1)) );
    for ( int i = 0;i<b.size();i+=nx) {
        b(i)=0;
    }

        vec a =1/ q(0,0,size(ng-(nx+1),1));
        for ( int j = nx-1;j<a.size();j+=nx) {
            a(j)=0;
        }

    s.diag(nx+1) = d ;
    s.diag(nx-1)=c;
    s.diag(-(nx-1))= b;
    s.diag(-(nx+1))  =a;

    return s/4;

}

sp_mat dxdy(vec &q, vec &r, double dx, double dy, int nx, int ny) {

}


