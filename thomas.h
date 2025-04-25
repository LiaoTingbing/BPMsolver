#pragma once

#include "common.h"

cx_vec thomas_algorithm(
    const cx_vec& a, const cx_vec& b, const cx_vec& c,
    const cx_vec& d);

cx_vec thomas_algorithm(
    const cx_vec& a, const cx_vec& b, const cx_vec& c,
    int p, const cx_vec& d);


 