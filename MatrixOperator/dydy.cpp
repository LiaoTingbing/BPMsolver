//
// Created by ltb on 25-4-19.
//

#include "dydy.h"

sp_mat dydy(double dy, int nx, int ny) {

    const int ng  = nx * ny ;
    sp_mat S(ng , ng ) ;

    vec a(ng -nx , fill::value(1));
    vec b(ng , fill::value(-2 ));

    S.diag(nx)=a;
    S.diag(-nx)=a;
    S.diag(0) =b;

    return S / dy / dy ;
}

sp_mat dydy(vec &q, double dy, int nx, int ny) {

    const int ng  = nx * ny ;
    sp_mat S(ng , ng ) ;

    vec z(nx);

    vec qr = q(0,0,size(ng-nx,1));
    vec ql = q(nx,0,size(ng-nx,1));
    vec qa = 1.0 / ( qr + ql  );

    vec qb = 1/ ( q + join_cols(z,qr)) + 1 /(q +  join_cols(ql,z));

    S.diag(nx)  = qa ;
    S.diag(-nx)  = qa;
    S.diag(0) = -qb;

    return S* 2 / dy / dy ;
}
