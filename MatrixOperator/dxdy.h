//
// Created by ltb on 25-4-19.
//

#ifndef DXDY_H
#define DXDY_H

#endif //DXDY_H

#pragma once

#include "../include/commom.h"


sp_mat dxdybase(vec & q , int nx ,int ny);

sp_mat dxdy( vec & q , vec & r ,double dx , double dy , int nx , int ny );