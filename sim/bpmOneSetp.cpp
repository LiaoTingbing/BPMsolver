//
// Created by ltb on 25-4-20.
//

#include "bpmOneSetp.h"


void onestepViolence(vec &ERX, vec &ERY, vec &ERZ, vec &ERXY, vec &ERYX, vec &ERX2, vec &ERY2, vec &ERZ2, vec &ERXY2,
                     vec &ERYX2, cx_vec &Exin, cx_vec &Eyin, cx_vec Et, double dx, double dy, double dz, double k0,
                     double n0,
                     double alpha, int nx, int ny,
                     int nz) {
        sp_mat P1 = calculateP(
                ERX, ERY, ERZ, ERXY, ERYX,
                dx, dy, k0, n0,
                nx, ny);

        sp_mat P2 = calculateP(
                ERX, ERY, ERZ, ERXY, ERYX,
                dx, dy, k0, n0,
                nx, ny);

    sp_cx_mat A = spdiags( 2.0*iu*n0*k0/dz * ones(2*nx*ny,1),ivec{0},2*nx*ny,2*nx*ny)
    -alpha*P2;

    sp_cx_mat B = (spdiags( 2.0*iu*n0*k0/dz*ones(2*nx*ny,1),ivec{0},2*nx*ny,2*nx*ny)
    +(1-alpha)*P1);

    cx_vec bb = B*Et;
    // cx_vec Eout = spsolve(A , B * Et);








}

sp_mat calculateP(
        vec &ERX, vec &ERY, vec &ERZ, vec &ERXY, vec &ERYX,
        double dx, double dy, double k0, double n0,
        int nx, int ny) {
        int ng = nx * ny;
        vec Iv(ng, fill::ones);
        sp_mat P(ng * 2, ng * 2);

        // Pxx
        P(0, 0, size(ng, ng)) =
                        dxdx(Iv, ERZ, ERX, nx, ny, dx, dy)
                        + dydy(Iv, Iv, Iv, nx, ny, dx, dy)
                        + spdiags(k0 * k0 * (ERX / eps0 - n0 * n0), ivec{0}, ng, ng);

        // Pyy
        P(ng, ng, size(ng, ng)) =
                        dxdx(Iv, Iv, Iv, nx, ny, dx, dy)
                        + dydy(Iv, ERZ, ERY, nx, ny, dx, dy)
                        + spdiags(k0 * k0 * (ERY / eps0 - n0 * n0), ivec{0}, ng, ng);

        // PXY
        P(0, ng, size(ng, ng)) =
                        dxdy(Iv, ERZ, ERY, nx, ny, dx, dy)
                        - dxdy(Iv, Iv, Iv, nx, ny, dx, dy)
                        + spdiags(ERXY / eps0 * k0 * k0, ivec{0}, ng, ng);

        // PYX
        P(ng, 0, size(ng, ng)) =
                        dydx(Iv, ERZ, ERX, nx, ny, dx, dy)
                        - dydx(Iv, Iv, Iv, nx, ny, dx, dy)
                        + spdiags(ERYX / eps0 * k0 * k0, ivec{0}, ng, ng);
        return P;
}

