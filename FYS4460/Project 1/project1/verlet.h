/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  10.02.2015
 *
 *  c++11 compiler
 */

#ifndef VERLET_H
#define VERLET_H

#include "acceleration.h"
#include "array.h"

template<size_t N, typename X, typename V, typename A, typename D> class VelocityVerlet;
template<typename X, typename V, typename A, typename D> class VelocityVerlet<2,X,V,A,D> : Acceleration<void,X,V,Array<void,A>,D> {
    protected: X* x;
    protected: V* v;
    protected: Array<void,A> a;
    public: inline VelocityVerlet(X& x, V& v, const D& f) : Acceleration<void,X,V,Array<void,A>,D>(f), x(&x), v(&v), a(x) {
        this->Run(x,v,a);           // Initialize
    }
    public: template<typename T> inline void Solve(T& dt) {
        T dt_div_2 = dt / 2;
        size_t n = a;
        a = static_cast<size_t>(*x);
        if(n != a)
            this->Run(*x,*v,a);     // Initialize
        for(size_t i = 0; i < *v; i++) {
            for(size_t j = 0; j < (*v)[i]; j++) {
                (*v)[i][j] += a[i][j] * dt_div_2;
                (*x)[i][j] += (*v)[i][j] * dt;
            }
        }
        this->Run(*x,*v,a);
        for(size_t i = 0; i < *v; i++) {
            for(size_t j = 0; j < (*v)[i]; j++)
                (*v)[i][j] += a[i][j] * dt_div_2;
        }
    }
    public: inline size_t& operator=(const size_t& len) {
        return a = len;
    }
    public: inline operator size_t&() {
        return a;
    }
    public: inline operator size_t() const {
        return a;
    }
};
template<typename X, typename V, typename A> class VelocityVerlet<2,X,V,A,void> : public VelocityVerlet<2,X,V,A,Delegate<void,void,X&,V&,Array<void,A>&>> {
    public: inline VelocityVerlet(X& x, V& v) :
        VelocityVerlet<2,X,V,A,Delegate<void,void,X&,V&,Array<void,A>&>>(x, v, delegate(&Acceleration<void,X,V,Array<void,A>,void>::Run)) {
    }
    public: template<typename T> inline void Solve(T& dt) {
        VelocityVerlet<2,X,V,A,Delegate<void,void,X&,V&,Array<void,A>&>>::Solve(dt);
    }
};
#endif // VERLET_H
