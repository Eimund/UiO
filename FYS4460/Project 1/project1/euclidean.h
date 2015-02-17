/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  17.02.2015
 *
 *  c++11 compiler
 */
#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H

#include "unit.h"
#include "vector.h"

struct Euclidean {
    template<typename X> static auto Length(const X& x) -> decltype(x[0])  {
        decltype(x[0]*x[0]) square = 0;
        for(size_t i = 0; i < x; i++)
            square += x[i]*x[i];
        return sqrt(square);
    }
    template<typename R, typename X, size_t D> static auto UnitVector(const R& r, const Vector<X,D>& x) -> Vector<decltype(x[0]/r),D> {
        Vector<decltype(x[0]/r),D> v;
        for(size_t i = 0; i < D; i++)
            v[i] = x[i]/r;
        return v;
    }
};

#endif // EUCLIDEAN_H
