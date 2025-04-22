#include "dydxbase.h"


sp_mat dydxbase(vec& q, int nx, int ny, double dx, double dy) {
    int ng = nx * ny;
    sp_mat s(ng, ng);

    vec a(ng), b(ng), c(ng), d(ng);
    vec iq = 1 / q;;

    int point3, point7;

    for (int i = 0; i < ng; i++) {
        point3 = i + nx;
        point7 = i - nx;


        if (i == 0) {
            // ◊Ûœ¬Ω«
            d(i) = iq(point3);
            // c(i) = iq(point3);
            // b(i) = iq(point7);
            // a(i) = iq(point7);
        }
        else if (i == nx - 1) {
            // ”“œ¬
            // d(i) = iq(point3);
            c(i) = iq(point3);
            // b(i) = iq(point7);
            // a(i) = iq(point7);
        }
        else if (i == ng - nx) {
            // ◊Û…œ
            // d(i) = iq(point3);
            // c(i) = iq(point3);
            b(i) = iq(point7);
            // a(i) = iq(point7);
        }
        else if (i == ng - 1) {
            // ”“…œ
            // d(i) = iq(point3);
            // c(i) = iq(point3);
            // b(i) = iq(point7);
            a(i) = iq(point7);
        }
        else if (i < nx) {
            // œ¬±ﬂ
            d(i) = iq(point3);
            c(i) = iq(point3);
            // b(i) = iq(point7);
            // a(i) = iq(point7);
        }
        else if (i % nx == 0) {
            //◊Û±ﬂ
            d(i) = iq(point3);
            // c(i) = iq(point3);
            b(i) = iq(point7);
            // a(i) = iq(point7);
        }
        else if (i % nx == nx - 1) {
            // ”“±ﬂ
            // d(i) = iq(point3);
            c(i) = iq(point3);
            // b(i) = iq(point7);
            a(i) = iq(point7);
        }
        else if (i > ng - nx) {
            // …œ±ﬂ
            // d(i) = iq(point3);
            // c(i) = iq(point3);
            b(i) = iq(point7);
            a(i) = iq(point7);
        }
        else {
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
