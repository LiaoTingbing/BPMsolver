
// #pragma once

#include "include/commom.h"
#include "include/Thomas.h"
#include "MatrixOperator/dxdx.h"
#include "MatrixOperator/dydy.h"
#include "MatrixOperator/dxdy.h"
#include "MatrixOperator/dydx.h"


int main() {


    double dx=0.004;
    double dy=0.004;

    int nx=4;
    int ny = 4;

    vec q(nx*ny , fill::ones);
    vec r(nx*ny , fill::ones);

      dxdybase(q,nx,ny)  ;

    return 0;
}