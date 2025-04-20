//
// Created by ltb on 25-4-20.
//

#ifndef DXDX_H
#define DXDX_H

#include "dxdxbase.h"

sp_mat dxdx(vec &p,
            vec &q,
            vec &r,
            int nx, int ny,
            double dx, double dy);


#endif //DXDX_H


