//
// Created by ltb on 25-4-19.
//

#ifndef DXDX_H
#define DXDX_H

#endif //DXDX_H


#pragma once
#include "../include/commom.h"

// 系数 1 1
sp_mat dxdx(double dx , int nx , int ny);

// 系数 1/q r
sp_mat dxdx(vec q , vec r , double dx , int nx , int ny);
