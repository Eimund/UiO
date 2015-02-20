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
#include "acceleration.h"
#include "array.h"
#include "chain.h"
#include "min_image.h"
#include "property.h"
#include "ref.h"
#include "space.h"
#include "template.h"
#include "unit.h"
#include "vector.h"

using namespace std;

template<typename X, typename V, typename A, typename F, typename T, size_t D> class PBC : public MinImage<T,D>, public Acceleration<void,Chain<Ref<typename RemoveRef<decltype(X()[0])>::type>>,Chain<Ref<typename RemoveRef<decltype(V()[0])>::type>> ,Chain<Ref<typename RemoveRef<decltype(A()[0])>::type>>,F> {
    public: typedef Chain<Ref<typename RemoveRef<decltype(X()[0])>::type>> CX;
    public: typedef Chain<Ref<typename RemoveRef<decltype(V()[0])>::type>> CV;
    public: typedef Chain<Ref<typename RemoveRef<decltype(A()[0])>::type>> CA;
    private: class _origin_ {
        private: Vector<Property<PBC<X,V,A,F,T,D>,T>,D> origin;
        public: inline Vector<Property<PBC<X,V,A,F,T,D>,T>,D>& operator=(const Vector<Property<PBC<X,V,A,F,T,D>,T>,D>& origin) {
            return this->origin = origin;
        }
        public: inline operator Vector<Property<PBC<X,V,A,F,T,D>,T>,D>&() {
            return origin;
        }
        public: inline operator Vector<Property<PBC<X,V,A,F,T,D>,T>,D>() const {
            return origin;
        }
        public: inline Property<PBC<X,V,A,F,T,D>,T>& operator[](size_t i) {
            return origin[i];
        }
        public: inline Property<PBC<X,V,A,F,T,D>,T> operator[](size_t i) const {
            return origin[i];
        }
    };
    public: _origin_ origin;
    private: class _range_ {
        private: Vector<Property<PBC<X,V,A,F,T,D>,T>,D> range;
        public: inline Vector<Property<PBC<X,V,A,F,T,D>,T>,D>& operator=(const Vector<Property<PBC<X,V,A,F,T,D>,T>,D>& range) {
            return this->range = range;
        }
        public: inline operator Vector<Property<PBC<X,V,A,F,T,D>,T>,D>&() {
            return range;
        }
        public: inline operator Vector<Property<PBC<X,V,A,F,T,D>,T>,D>() const {
            return range;
        }
        public: inline Property<PBC<X,V,A,F,T,D>,T>& operator[](size_t i) {
            return range[i];
        }
        public: inline Property<PBC<X,V,A,F,T,D>,T> operator[](size_t i) const {
            return range[i];
        }
    };
    public: _range_ range;
    private: Vector<Array<PBC<X,V,A,F,T,D>,T>,D> N;
    private: Vector<Property<PBC<X,V,A,F,T,D>,size_t>,D> n;
    private: template<typename O, size_t S> struct SetSize : SetSize<O,S-1> {
        Delegate<O,void,T*,size_t,size_t> func1;
        Delegate<O,void,T&,const T&> func2;
        Delegate<O,void,size_t&,const size_t&> func3;
        SetSize(O* owner) : SetSize<O,S-1>(owner), func1(delegate(*owner, &O::template ImageWindows<S-1>)),
            func2(delegate(*owner, &O::template SetRange<S-1>)), func3(delegate(*owner, &O::template SetN<S-1>)) {
            owner->N[S-1].init = &func1;
            owner->N[S-1] = 2;
            owner->origin[S-1].func = &func2;
            owner->range[S-1].func = &func2;
            owner->n[S-1].func = &func3;
        }
    };
    private: template<typename O> struct SetSize<O,0> {
        O* owner;
        SetSize(O *owner) : owner(owner) {
        }
    };
    private: SetSize<PBC<X,V,A,F,T,D>,D> size;
    private: template<typename Z, size_t S> struct Within {
        template<typename CX, typename CV, typename CA, typename DX, typename DV, typename DA> static inline void Check(const Z& owner, CX& cx, CV& cv, CA& ca, const DX& x, const DV& v, const DA& a) {
            size_t n1 = owner.N[S-1]-1;
            for(size_t i = 0; i < n1; i++) {
                if(owner.N[S-1][i] <= x[S-1] && owner.N[S-1][i+1] > x[S-1]) {
                    Within<Z,S-1>::Check(owner, cx[i], cv[i], ca[i], x, v, a);
                    break;
                }
            }
        }
        template<typename CX, typename CV, typename CA> static inline void Run(Z& owner, CX& x, CV& v, CA& a, size_t n[S]) {
            for(size_t i = 0; i < n[0]; i++)
                Within<Z,S-1>::Run(owner, x[i], v[i], a[i], &n[1]);
        }
    };
    private: template<typename Z> struct Within<Z,0> {
        template<typename CX, typename CV, typename CA, typename DX, typename DV, typename DA> static inline void Check(const Z&, CX& cx, CV& cv, CA& ca, const DX& x, const DV& v, const DA& a) {
            cx += x;
            cv += v;
            ca += a;
        }
        template<typename CX, typename CV, typename CA> static inline void Run(Z& owner, CX& x, CV& v, CA& a, size_t[0]) {
            owner.Run(x,v,a);
        }
    };
    public: PBC(F& acc) : Acceleration<void,CX,CV,CA,F>(acc), size(this) {
    }
    public: inline Property<PBC<X,V,A,F,T,D>,size_t>& operator[](size_t i) {
        return n[i];
    }
    public: template<typename U> void Boundary(U& pos) {
        Vector<T,D> dx;
        T buf, buf2;
        for(size_t i = 0; i < D; i++)
            dx[i] = range[i] - origin[i];
        for(size_t i = 0; i < pos; i++) {
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
    private: template<size_t S> void ImageWindows(T* arr, size_t, size_t n) {
        auto dx = range[S]-origin[S];
        n--;
        this->L[S] = dx;// /n;
        for(size_t i = 0; i <= n; i++)
            arr[i] = origin[S] + dx * (static_cast<typename RemoveTemplate<T>::type>(i) / n);
    }
    private: template<size_t S> void SetRange(T& old, const T& ny) {
        old = ny;
        size_t n = N[S];
        N[S] = 0;
        N[S] = n;
    }
    private: template<size_t S> void SetN(size_t& old, const size_t& ny) {
        if(ny == 0)
            old = ny+2;
        else
            old = ny+1;
        N[S] = 0;
        N[S] = old;
    }
    public: void NeighbourList(X& x, V& v, A& a) {
        size_t len[D];
        for(size_t i = 0; i < D; i++)
            len[i] = N[i]-1;
        auto arr_x = Space<CX,D>::Allocate(len);
        auto arr_v = Space<CV,D>::Allocate(len);
        auto arr_a = Space<CA,D>::Allocate(len);

        for(size_t i = 0; i < x; i++) { // Initialize forces to zero
            for(size_t j = 0; j < D; j++)
                a[i][j] = 0;
        }

        for(size_t i = 0; i < x; i++)   // Find the correct box for each particle
            Within<PBC<X,V,A,F,T,D>,D>::Check(*this, arr_x, arr_v, arr_a, x[i], v[i], a[i]);

        if(len[0] > 1) {
            auto x0 = arr_x[0][0][0];
            auto x1 = arr_x[0][0][1];
            auto x2 = arr_x[0][1][0];
            auto x3 = arr_x[0][1][1];
            auto x4 = arr_x[1][0][0];
            auto x5 = arr_x[1][0][1];
            auto x6 = arr_x[1][1][0];
            auto x7 = arr_x[1][1][1];
            auto x8 = 0;
        }

        Within<PBC<X,V,A,F,T,D>,D>::Run(*this, arr_x, arr_v, arr_a, len);  // Run simulation on particles

        Space<CA,D>::Deallocate(arr_a, len);
        Space<CV,D>::Deallocate(arr_v, len);
        Space<CX,D>::Deallocate(arr_x, len);
    }
};

template<typename X, typename V, typename A, typename T, size_t D> class PBC<X,V,A,void,T,D> {
    public: typedef Chain<Ref<typename RemoveRef<decltype(X()[0])>::type>> CX;
    public: typedef Chain<Ref<typename RemoveRef<decltype(V()[0])>::type>> CV;
    public: typedef Chain<Ref<typename RemoveRef<decltype(A()[0])>::type>> CA;
};

#endif // PBC_H
