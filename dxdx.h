#pragma once

#include "common.h"
#include "dxdxbase.h"

sp_mat dxdx(vec& p,
    vec& q,
    vec& r,
    int nx, int ny,
    double dx, double dy);
