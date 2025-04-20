//
// Created by ltb on 25-4-19.
//

#include "dxdxbase.h"

#include <ranges>


sp_mat dxdxbase(vec &q, int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat s(ng, ng);

    vec a(ng), b(ng), c(ng);

    int point1, point2, point3;
    for (int i = 0; i < ng; i++) {
        point2 = i;
        point1 = i - 1;
        point3 = i + 1;

        if (i % nx == 0) {
            // 左边界
            // a(i) = 1 / ( q(point1) + q(point2) ) ;
            b(i) = -(1 / (0 + q(point2)) + 1 / (q(point2) + q(point3)));
            c(i) = 1 / (q(point2) + q(point3));
        } else if (i % nx == nx - 1) {
            // 右边界
            a(i) = 1 / (q(point1) + q(point2));
            b(i) = -(1 / (q(point1) + q(point2)) + 1 / (q(point2) + 0));
            // c(i) = 1 / (q(point2) + q(point3));
        } else {
            a(i) = 1 / (q(point1) + q(point2));
            b(i) = -(1 / (q(point1) + q(point2)) + 1 / (q(point2) + q(point3)));
            c(i) = 1 / (q(point2) + q(point3));
        }
    }

    s.diag(-1) = a(1, 0, size(ng - 1, 1));
    s.diag(0) = b;
    s.diag(1) = c(0, 0, size(ng - 1, 1));

    s = 2 * s / dx / dx;
    return s;
}
