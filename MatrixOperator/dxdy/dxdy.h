//
// Created by ltb on 25-4-20.
//

#ifndef DXDY_H
#define DXDY_H


#endif //DXDY_H

#pragma once
#include "dxdybase.h"

inline sp_mat dxdy(vec &p,
                   vec &q,
                   vec &r,
                   int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat P(ng, ng), R(ng, ng);

    P.diag(0) = p;
    R.diag(0) = r;
    return P * dxdybase(q, nx, ny, dx, dy) * R;
}
