/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  14.02.2015
 *
 *  c++11 compiler
 */
#ifndef LJPOTENTIAL_H
#define LJPOTENTIAL_H

#include "delegate.h"
#include "unit.h"

template<size_t D,typename X, typename A, typename M, typename F, typename E, typename C1, typename C2,typename C3> class LJ_Potential {
    private: X sigma;
    private: M mass;
    private: E epsilon;
    private: F force;
    private: A acc;
    public: Delegate<C1,Vector<X,D>,const Vector<X,D>&,const Vector<X,D>&>* dist_vec;
    public: Delegate<C2,X,const Vector<X,D>&>* dist;
    public: Delegate<C3,Vector<decltype(sigma/sigma),D>,const X&, const Vector<X,D>&>* unit_vec;
    public: inline LJ_Potential(X sigma, M mass, E epsilon) : sigma(sigma), mass(mass), epsilon(epsilon),
        force(4*epsilon/sigma), acc(force/mass), dist_vec(nullptr), dist(nullptr), unit_vec(nullptr) {
    }
    public: inline X& operator=(const X& sigma) {
        this->sigma = sigma;
        Update();
        return this->sigma;
    }
    public: inline M& operator=(const M& mass) {
        this->mass = mass;
        Update();
        return this->mass;
    }
    public: inline E& operator=(const E& epsilon) {
        this->epsilon = epsilon;
        Update();
        return this->epsilon;
    }
    public: inline A Acceleration(const X& r) {
        auto x = sigma/r;
        auto p = x*x*x;
        p *= p;
        return acc*x*(p-p*p);
    }
    public: template<typename X1, typename A1> inline void Acceleration(X1& x, A1& a) {
        for(size_t i = 0; i < x; i++) {
            for(size_t j = 0; j < D; j++)
                a[i][j] = 0;
        }
        for(size_t i = 0; i < x; i++) {
            for(size_t j = i+1; j < x; j++) {
                auto r_vec = (*dist_vec)(x[i], x[j]);
                auto r = (*dist)(r_vec);
                if(r != 0) {
                    auto u_vec = (*unit_vec)(r, r_vec);
                    auto aks = Acceleration(r);
                    for(size_t k = 0; k < D; k++) {
                        auto ak = aks*u_vec[k];
                        a[i][k] -= ak;
                        a[j][k] += ak;
                    }
                }
            }
        }
    }
    private: inline void Update() {
        force = 4*epsilon/sigma;
        acc = force/mass;
    }

};


#endif // LJPOTENTIAL_H
