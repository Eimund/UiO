/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  13.02.2015
 *
 *  c++11 compiler
 */

#ifndef PBC_H
#define PBC_H

#include <cmath>
#include "unit.h"
#include "vector.h"

using namespace std;

template<typename T, size_t D> struct PBC {
    Vector<T,D> origin;
    Vector<T,D> range;
    template<typename U> void Boundary(U& pos) {
        Vector<T,D> dx;
        T buf;
        for(size_t i = 0; i < D; i++)
            dx[i] = range[i] - origin[i];
        for(size_t i = 0; i < static_cast<size_t>(pos); i++) {
            for(size_t j = 0; j < D; j++) {
                if(pos[i][j] < origin[j]) {
                    buf = fmod(dx[j],origin[j]-pos[i][j]);
                    pos[i][j] = range[j] - buf;
                } else if(pos[i][j] > range[j]) {
                    buf = fmod(dx[j],pos[i][j]-range[j]);
                    pos[i][j] = origin[j] + buf;
                }
            }
        }
    }
};

#endif // PBC_H
