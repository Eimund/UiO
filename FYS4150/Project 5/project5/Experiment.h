/*
 *  FYS4150 - Computational Physics - Project 5
 *
 *  Written by:         Eimund Smestad
 *  Candiate number:    68
 *
 *  30.11.2014
 *
 *  c++11 compiler
 */

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "Array.h"
#include "Delegate.h"

template<typename V, typename T, unsigned int D, typename C, typename... P> typename Pointer<T,D>::Type Experiment(V N, T* x[D], unsigned int n[D], Delegate<C,Chain<Vector<T,D>>,P...> experiment, P... arg) {

    auto space = Space<V,D>::Allocate(n);

    for(int i = 0; i < N; i++) {
        Chain<Vector<T,D>> states = experiment(arg...);
        Space<T,D>::template Map<V>(space, states, x, n);
    }

    auto space2 = Space<T,D>::template Normalize<V>(space, n, Space<V,D>::Max(space, n));
    Space<V,D>::Deallocate(space, n);
    return space2;
}

template<typename V, typename T, unsigned int D, typename C, typename... P> typename Pointer<T,D>::Type Experiment2(V N, T* x[D], unsigned int n[D], Delegate<C,void,typename Pointer<T,D>::Type,T*[D],unsigned int[D],P...> experiment, P... arg) {

    auto s1 = Space<T,D>::Allocate(n);
    auto s2 = Space<T,D>::Allocate(n);

    for(int i = 0; i < N; i++) {
        experiment(s1, x, n, arg...);
        Space<T,D>::Add(s1, s2, n);
    }

    if(N > 1)
        Space<T,D>::Normalize(s1, n, static_cast<T>(N));
    Space<T,D>::Deallocate(s2, n);
    return s1;
}

template<typename T> void Metropolis(T T1, T& w1, T T2, T& w2) {
    T a1 = T1 * w1;
    T a2 = T2 * w2;
    T A;
    if(a1 <= a2)
        A = 1;
    else
        A = a2/a1;
    w1 += A*a2 - A*a1;
    w2 += A*a1 - A*a2;
}

#endif // EXPERIMENT_H
