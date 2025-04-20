//
// Created by ltb on 25-4-19.
//

#include "dydxbase.h"

sp_mat dydxbase(vec &q, int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat s(ng, ng);

    vec a(ng), b(ng), c(ng), d(ng);
    vec iq = 1 / q;;

    int point3, point7;

    for (int i = 0; i < ng; i++) {
        point3 = i + nx;
        point7 = i - nx;


        if (i == 0) {
            // 左下角
            d(i) = iq(point3);
            // c(i) = iq(point3);
            // b(i) = iq(point7);
            // a(i) = iq(point7);
        } else if (i == nx - 1) {
            // 右下
            // d(i) = iq(point3);
            c(i) = iq(point3);
            // b(i) = iq(point7);
            // a(i) = iq(point7);
        } else if (i == ng - nx) {
            // 左上
            // d(i) = iq(point3);
            // c(i) = iq(point3);
            b(i) = iq(point7);
            // a(i) = iq(point7);
        } else if (i == ng - 1) {
            // 右上
            // d(i) = iq(point3);
            // c(i) = iq(point3);
            // b(i) = iq(point7);
            a(i) = iq(point7);
        } else if (i < nx) {
            // 下边
            d(i) = iq(point3);
            c(i) = iq(point3);
            // b(i) = iq(point7);
            // a(i) = iq(point7);
        } else if (i % nx == 0) {
            //左边
            d(i) = iq(point3);
            // c(i) = iq(point3);
            b(i) = iq(point7);
            // a(i) = iq(point7);
        } else if (i % nx == nx - 1) {
            // 右边
            // d(i) = iq(point3);
            c(i) = iq(point3);
            // b(i) = iq(point7);
            a(i) = iq(point7);
        } else if (i > ng - nx) {
            // 上边
            // d(i) = iq(point3);
            // c(i) = iq(point3);
            b(i) = iq(point7);
            a(i) = iq(point7);
        } else {
            d(i) = iq(point3);
            c(i) = iq(point3);
            b(i) = iq(point7);
            a(i) = iq(point7);
        }
    }

    s.diag(nx + 1) = d(0, 0, size(ng - (nx + 1), 1));
    s.diag(nx - 1) = -c(0, 0, size(ng - (nx - 1), 1));
    s.diag(-(nx - 1)) = -b(nx - 1, 0, size(ng - (nx - 1), 1));
    s.diag(-(nx + 1)) = a(nx + 1, 0, size(ng - (nx + 1), 1));
    s = s / 4 / dx / dy;
    return s;
}
