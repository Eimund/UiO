/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Array.h"
#include "Delegate.h"
#include "Property.h"

template<typename T, unsigned int D> class Boundary {
    private: typedef Boundary<T, D> THIS;
    protected: T* x[D+1];
    public: Property<PropertyType::ReadOnly, THIS, PropertySet<THIS, T>, T> lower[D+1];
    public: Property<PropertyType::ReadOnly, THIS, PropertySet<THIS, T>, T> upper[D+1];
    public: Property<PropertyType::ReadOnly, THIS, ArrayLength<THIS, T>, int> n[D+1];
    public: Boundary(T lower[D+1], T upper[D+1], int n[D+1]) {
        Init<0>(lower[0], upper[0], n[0]);
        Init<1>(lower[1], upper[1], n[1]);
    }
    private: template<unsigned int DIM> void Init(T lower, T upper, int n) {
        this->lower[DIM] = Property<PropertyType::ReadOnly, THIS, PropertySet<THIS, T>, T>(PropertySet<THIS, T>(lower, Delegate<THIS, T,T>(this, &THIS::Lower<DIM>)));
        this->upper[DIM] = Property<PropertyType::ReadOnly, THIS, PropertySet<THIS, T>, T>(PropertySet<THIS, T>(upper, Delegate<THIS, T,T>(this, &THIS::Upper<DIM>)));
        this->n[DIM] = Property<PropertyType::ReadOnly, THIS, ArrayLength<THIS, T>, int>(ArrayLength<THIS, T>(x[DIM], n, Delegate<THIS, void, T*, unsigned int, unsigned int>(this, &THIS::N<DIM>)));
    }
    private: template<unsigned int DIM> T Lower(T lower) {
        T step = (upper[DIM]-lower)/(n[DIM]-1);
        for(int j = 0; j < n[DIM]; j++)
            x[DIM][j] = lower + step*j;
        return lower;
    }
    private: template<unsigned int DIM> void N(T* array, unsigned int, unsigned n) {
        T step = (upper[DIM]-lower[DIM])/(n-1);
        for(int j = 0; j < n; j++)
            array[j] = lower[DIM] + step*j;
    }
    private: template<unsigned int DIM> T Upper(T upper) {
        T step = (upper-lower[DIM])/(n[DIM]-1);
        for(int j = 0; j < n[DIM]; j++)
            x[DIM][j] = lower[DIM] + step*j;
        return upper;
    }
};

#endif // BOUNDARY_H
