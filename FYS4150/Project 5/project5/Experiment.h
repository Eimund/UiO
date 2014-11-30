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

#endif // EXPERIMENT_H
