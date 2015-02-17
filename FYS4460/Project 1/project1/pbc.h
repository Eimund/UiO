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

template<typename X, typename V, typename A, typename F, typename T, size_t D> class PBC : public MinImage<T,D>, public Acceleration<X,V,A,F> {
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
    SetSize<PBC<X,V,A,F,T,D>,D> size;
    private: template<typename U, typename W, typename Y, size_t S> struct Within {
        static void Check(const U& chain, const W& data, const Y& N) {
            for(size_t i = 1; i < N[S-1]; i++) {
                if(N[S-1][i-1] >= data[S-1] && N[S-1][i] <= data[S-1])
                    Within<U,W,Y,S-1>::Check(chain[i], data, N);
            }
        }
    };
    private: template<typename U, typename W, typename Y> struct Within<U,W,Y,0> {
        static void Check(const U& chain, const V& data, const W& N) {

        }
    };
    public: PBC(F& acc) : Acceleration<X,V,A,F>(acc), size(this) {
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
        this->L[S] = dx/n;
        for(size_t i = 0; i <= n; i++)
            arr[i] = origin[S] + dx * (static_cast<typename RemoveTemplate<T>::type>(i) / n);
    }
    public: bool Inside(const Vector<T,D>& x1, const Vector<T,D>& x2) {
        return false;
    }
    private: template<size_t S> void SetRange(T& old, const T& ny) {
        old = ny;
        size_t n = N[S];
        N[S] = 0;
        N[S] = n;
    }
    private: template<size_t S> void SetN(size_t& old, const size_t& ny) {
        old = ny+1;
        N[S] = 0;
        N[S] = old;
    }
    public: inline Property<PBC<X,V,A,F,T,D>,size_t>& operator[](size_t i) {
        return n[i];
    }
    public: void NeighbourList(X& x, V& v, A& a) {
        size_t len[D];
        for(size_t i = 0; i < D; i++)
            len[i] = N[i]-1;
        auto arr_x = Space<Chain<Ref<typename RemoveRef<decltype(x[0][0])>::type>>,D>::Allocate(len);
        auto arr_v = Space<Chain<Ref<typename RemoveRef<decltype(x[0][0])>::type>>,D>::Allocate(len);
        auto arr_a = Space<Chain<Ref<typename RemoveRef<decltype(x[0][0])>::type>>,D>::Allocate(len);

        this->Run(x,v,a);

        Space<Chain<Ref<typename RemoveRef<decltype(x[0][0])>::type>>,D>::Deallocate(arr_a, len);
        Space<Chain<Ref<typename RemoveRef<decltype(x[0][0])>::type>>,D>::Deallocate(arr_v, len);
        Space<Chain<Ref<typename RemoveRef<decltype(x[0][0])>::type>>,D>::Deallocate(arr_x, len);
    }
};

#endif // PBC_H
