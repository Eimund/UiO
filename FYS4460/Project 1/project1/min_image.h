/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  16.02.2015
 *
 *  c++11 compiler
 */
#ifndef MIN_IMAGE_H
#define MIN_IMAGE_H

#include "unit.h"
#include "vector.h"

template<typename T, size_t D> class MinImage {
    protected: Vector<T,D> L{0};
    public: template<typename X1, typename X2> Vector<T,D> Distance(const Vector<X1,D>& x1, const Vector<X2,D>& x2) {
        Vector<T,D> r(0);
        T buf1, buf2;
        for(size_t i = 0; i < D; i++) {
            r[i] = x1[i]-x2[i];
            buf1 = r[i] + L[i];
            buf2 = r[i] - L[i];
            if(abs(buf1) < abs(r[i]))
                r[i] = buf1;
            if(abs(buf2) < abs(r[i]))
                r[i] = buf2;
        }
        return r;
    }
};

#endif // MIN_IMAGE_H
