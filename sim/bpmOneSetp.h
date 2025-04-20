//
// Created by ltb on 25-4-20.
//

#ifndef BPMONESETP_H
#define BPMONESETP_H

#include "../MatrixOperator/dxdx/dxdx.h"
#include "../MatrixOperator/dxdy/dxdy.h"
#include "../MatrixOperator/dydx/dydx.h"
#include "../MatrixOperator/dydy/dydy.h"



#include "../function/commom.h"


void onestepViolence(
    vec &ERX, vec &ERY, vec &ERZ, vec &ERXY, vec & ERYX,
    vec &ERX2, vec &ERY2, vec &ERZ2, vec &ERXY2, vec &  ERYX2,
    cx_vec &Exin, cx_vec & Eyin, cx_vec Et,
    double dx, double dy, double dz, double k0, double n0,double alpha,
    int nx, int ny, int nz
);


sp_mat calculateP(
    vec &ERX, vec &ERY, vec &ERZ, vec &ERXY, vec & ERYX,
    double dx, double dy, double k0, double n0,
    int nx, int ny
) ;

#endif //BPMONESETP_H