#pragma once


#include "dydybase.h"


inline sp_mat dydy(vec& p,
    vec& q,
    vec& r,
    int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat P(ng, ng), R(ng, ng);

    P.diag(0) = p;
    R.diag(0) = r;
    return P * dydybase(q, nx, ny, dx, dy) * R;
}
