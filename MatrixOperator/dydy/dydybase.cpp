//
// Created by ltb on 25-4-19.
//

#include "dydybase.h"

sp_mat dydybase(vec &q, int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat s(ng, ng);
    vec a(ng), b(ng), c(ng);

    int point1, point2, point3;

    for (int i = 0; i < ng; i++) {
        point2 = i;
        point1 = i - nx;
        point3 = i + nx;

        if (i < nx) {
            // 下边界
            // a(i) = 1 / ( q(point1) + q(point2) ) ;
            b(i) = -(1 / (0 + q(point2) )+ 1 / (q(point2) + q(point3)));
            c(i) = 1 / (q(point2) + q(point3));
        } else if (i > ng - nx - 1) {
            // 上边界
            a(i) = 1 / (q(point1) + q(point2));
            b(i) = -(1 / (q(point1) + q(point2) )+ 1 / (q(point2) + 0));

            // c(i) = 1 / (q(point2) + q(point3));
        } else {
            a(i) = 1 / (q(point1) + q(point2));
            b(i) = -(1 / (q(point1) + q(point2) )+ 1 / (q(point2) + q(point3)));
            c(i) = 1 / (q(point2) + q(point3));
        }
    }

    s.diag(nx) = c(0, 0, size(ng - nx, 1));
    s.diag(0) = b;
    s.diag(-nx) = a(nx, 0, size(ng - nx, 1));

    s = 2 * s / dy / dy;
    return s;
}
