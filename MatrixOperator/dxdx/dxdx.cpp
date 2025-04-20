//
// Created by ltb on 25-4-20.
//

#include "dxdx.h"


sp_mat dxdx(vec &p,
            vec &q,
            vec &r,
            int nx, int ny,
            double dx, double dy) {
    int ng = nx * ny;
    // cout<<ng;

    // sp_mat P(ng, ng), R(ng, ng);
    // P.diag(0) = p;
    // R.diag(0) = R;
    sp_mat P = spdiags(p, ivec{0}, ng, ng);
    sp_mat R = spdiags(r, ivec{0}, ng, ng);

    return P * dxdxbase(q, nx, ny, dx, dy) * R;
}
