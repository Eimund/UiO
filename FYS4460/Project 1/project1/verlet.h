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

#include "array.h"
#include "delegate.h"

template<typename X, typename V, typename A, typename D> class VelocityVerletFunction {

};

template<typename X, typename V, typename A> class VelocityVerletFunction<X,V,A,void> {
    public: static inline void Run(X&, V&, A&) {
    }
};
template<typename X, typename V, typename A, typename C> struct VelocityVerletFunction<X,V,A,Delegate<C,void,A&>> {
    private: Delegate<C,void,A&> f;
    protected: inline VelocityVerletFunction(const Delegate<C,void,A&>& f) : f(f) {
    }
    protected: inline void Run(X&, V&, A& a) {
        f(a);
    }
};
template<typename X, typename V, typename A, typename C> struct VelocityVerletFunction<X,V,A,Delegate<C,void,X&,A&>> {
    private: Delegate<C,void,X&,A&> f;
    protected: inline VelocityVerletFunction(const Delegate<C,void,X&,A&>& f) : f(f) {
    }
    protected: inline void Run(X& x, V&, A& a) {
        f(x,a);
    }
};
template<typename X,typename V, typename A, typename C> struct VelocityVerletFunction<X,V,A,Delegate<C,void,V&,A&>> {
    private: Delegate<C,void,V&,A&> f;
    protected: inline VelocityVerletFunction(const Delegate<C,void,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(X&, V& v, A& a) {
        f(v,a);
    }
};
template<typename X,typename V, typename A, typename C> struct VelocityVerletFunction<X,V,A,Delegate<C,void,X&,V&,A&>> {
    private: Delegate<C,void,X&,V&,A&> f;
    protected: inline VelocityVerletFunction(const Delegate<C,void,X&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(X& x, V& v, A& a) {
        f(x,v,a);
    }
};

template<size_t N, typename X, typename V, typename A, typename D> class VelocityVerlet;
template<typename X, typename V, typename A, typename D> class VelocityVerlet<2,X,V,A,D> : VelocityVerletFunction<X,V,Array<void,A>,D> {
    private: X* x;
    private: V* v;
    private: Array<void,A> a;
    public: inline VelocityVerlet(X& x, V& v, const D& f) : VelocityVerletFunction<X,V,Array<void,A>,D>(f), x(&x), v(&v), a(x) {
        f(x,v,a);
    }
    public: template<typename T> inline void Solve(T& dt) {
        T dt_div_2 = dt / 2;
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
};
template<typename X, typename V, typename A> class VelocityVerlet<2,X,V,A,void> : public VelocityVerlet<2,X,V,A,Delegate<void,void,X&,V&,Array<void,A>&>> {
    public: inline VelocityVerlet(X& x, V& v) :
        VelocityVerlet<2,X,V,A,Delegate<void,void,X&,V&,Array<void,A>&>>(x, v, delegate(&VelocityVerletFunction<X,V,Array<void,A>,void>::Run)) {
    }
    public: template<typename T> inline void Solve(T& dt) {
        VelocityVerlet<2,X,V,A,Delegate<void,void,X&,V&,Array<void,A>&>>::Solve(dt);
    }
};
#endif // VERLET_H
