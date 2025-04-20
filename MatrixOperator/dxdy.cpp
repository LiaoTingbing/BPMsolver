//
// Created by ltb on 25-4-19.
//

#include "dxdy.h"


sp_mat dxdybase(vec &q, int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat s(ng, ng);

    vec a(ng), b(ng), c(ng), d(ng);

    int point0, point2, point4, point6, point8;
    for (int i = 0; i < ng; i++) {
        // point0 = i ;
        // point1 = point0 +1 ;
        point2 = i + 1 + nx;
        // point3 = point0 + nx ;
        point4 = i - 1 + nx;
        // point5 = point0-1;
        point6 = i - 1 - nx;
        // point7 = point0 - nx ;
        point8 = i + 1 - nx;
        // cout<<i<<endl;

        vec iq = 1/q; ;

        if (i == 0) {
            // 左下角
            d(i) = iq(point2);
            // c(i) = -iq(point4);
            // b(i) = -iq(point8);
            // a(i) = iq(point6);
        } else if (i == nx - 1) {
            // 右下
            // d(i) = iq(point2);
            c(i) = -iq(point4);
            // b(i) = -iq(point8);
            // a(i) = iq(point6);
        } else if (i == ng - nx) {
            // 左上
            // d(i) = iq(point2);
            // c(i) = -iq(point4);
            b(i) = -iq(point8);
            // a(i) = iq(point6);
        } else if (i == ng - 1) {
            // 右上
            // d(i) = iq(point2);
            // c(i) = -iq(point4);
            // b(i) = -iq(point8);
            a(i) = iq(point6);
        } else if (i < nx) {
            // 下边
            d(i) = iq(point2);
            c(i) = -iq(point4);
            // b(i) = -iq(point8);
            // a(i) = iq(point6);
        } else if (i % nx == 0) {
            //左边
            d(i) = iq(point2);
            // c(i) = -iq(point4);
            b(i) = -iq(point8);
            // a(i) = iq(point6);
        } else if (i % nx == nx - 1) {
            // 右边
            // d(i) = iq(point2);
            c(i) = -iq(point4);
            // b(i) = -iq(point8);
            a(i) = iq(point6);
        } else if (i > ng - nx) {
            // 上边
            // d(i) = iq(point2);
            // c(i) = -iq(point4);
            b(i) = -iq(point8);
            a(i) = iq(point6);
        } else {
            d(i) = iq(point2);
            c(i) = -iq(point4);
            b(i) = -iq(point8);
            a(i) = iq(point6);
        }
    }

    s.diag(nx + 1) =   d(0, 0, size(ng - (nx + 1), 1));
    s.diag(nx - 1) =   c(0, 0, size(ng - (nx - 1), 1));
    s.diag(-(nx - 1)) =  b(nx - 1, 0, size(ng - (nx - 1), 1));
    s.diag(-(nx + 1)) =  a(nx + 1, 0, size(ng - (nx + 1), 1));
    s = s / 4 / dx / dy;
    return s;
}
