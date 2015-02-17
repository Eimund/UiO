/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  17.02.2015
 *
 *  c++11 compiler
 */
#ifndef SPACE_H
#define SPACE_H

#include <iostream>
#include <fstream>
#include "pointer.h"

using namespace std;

template<typename T, size_t D> struct Space {
    static void Add(typename Pointer<T,D>::Type s1, typename Pointer<T,D>::Type s2, size_t n[D]) {
        for(size_t i = 0; i < n[0]; i++)
            Space<T,D-1>::Add(s1[i], s2[i], &n[1]);
    }
    static inline typename Pointer<T,D>::Type Allocate(size_t n[D]) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(size_t i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::Allocate(&n[1]);
        return array;
    }
    static void ArrayFile(ofstream& file, typename Pointer<T,D>::Type array, size_t n[D]) {
        for(size_t i = 0; i < n[0]; i++) {
            if(D == 1 && i)
                file << '\t';
            else if(D > 1 && i)
                file << endl;
            Space<T,D-1>::ArrayFile(file, array[i], &n[1]);
        }
    }
    template<typename V> static typename Pointer<T,D>::Type Cast(typename Pointer<V,D>::Type space, size_t n[D]) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(size_t i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::Cast(space[i], &n[1]);
        return array;
    }
    static void Deallocate(typename Pointer<T,D>::Type array, size_t n[D]) {
        for(size_t i = 0; i < n[0]; i++)
            Space<T,D-1>::Deallocate(array[i], &n[1]);
        delete [] array;
    }
    static void DeRange(T** range) {
        for(size_t i = 0; i < D; i++)
            delete [] range[i];
        delete [] range;
    }
    template<typename V> static typename Pointer<V,D>::Type Map(Chain<Vector<T,D>>& chain, T* range[D], size_t n[D]) {
        auto space = Space<V,D>::Allocate(n);
        auto p = chain.owner;
        for(size_t i = 0; i < chain.N; i++) {
            p = p->next;
            Space<T,D>::Mapping<V>(space, p->e.e, range, n);
        }
        return space;
    }
    template<typename V> static void Map(typename Pointer<V,D>::Type& space, Chain<Vector<T,D>>& chain, T* range[D], size_t n[D]) {
        auto p = chain.owner;
        for(size_t i = 0; i < chain.N; i++) {
            p = p->next;
            Space<T,D>::Mapping<V>(space, p->e.e, range, n);
        }
    }
    template<typename V> static void Mapping(typename Pointer<V,D>::Type space, T element[D], T* range[D], size_t n[D]) {
        size_t m = n[0]-1;
        if(element[0] >= range[0][0] || element[0] <= range[0][m]) {
            for(size_t i = 0; i < m; i++) {
                if(element[0] < (range[0][i]+range[0][i+1])/2) {
                    Space<T,D-1>::template Mapping<V>(space[i], &element[1], &range[1], &n[1]);
                    return;
                }
            }
            Space<T,D-1>::template Mapping<V>(space[m], &element[1], &range[1], &n[1]);
        }
    }
    static T Max(typename Pointer<T,D>::Type space, size_t n[D]) {
        T val, max = 0;
        for(size_t i = 0; i < n[0]; i++) {
            val = Space<T,D-1>::Max(space[i], &n[1]);
            if(val > max)
                max = val;
        }
        return max;
    }
    static void Normalize(typename Pointer<T,D>::Type space, size_t n[D], T val) {
        for(size_t i = 0; i < n[0]; i++)
            Space<T,D-1>::Normalize(space[i], &n[1], val);
    }
    template<typename V> static typename Pointer<T,D>::Type Normalize(typename Pointer<V,D>::Type space, size_t n[D], V val) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(size_t i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::template Normalize<V>(space[i], &n[1], val);
        return array;
    }
    static T** Range(T lower[D], T upper[D], size_t n[D]) {
        T** array = new T*[D];
        for(size_t i = 0; i < D; i++) {
            T step = (upper[i]-lower[i])/(n[i]-1);
            array[i] = new T[n[i]];
            array[i][0] = lower[i];
            for(size_t j = 1; j < n[i]; j++)
                array[i][j] = array[i][j-1] + step;
        }
        return array;
    }
    static void RangeFile(ofstream& file, T* range[D], size_t n[D]) {
        for(size_t i = 0; i < D; i++) {
            if(i)
                file << endl;
            Space<T,1>::ArrayFile(file, range[i], &n[i]);
        }
    }
    static void ToFile(ofstream& file, T* range[D], typename Pointer<T,D>::Type array, size_t n[D]) {
        file << D << endl;
        RangeFile(file, range, n);
        file << endl;
        ArrayFile(file, array, n);
    }
};
template<typename T> struct Space<T,0> {
    static void Add(T& s1, T& s2, size_t[0]) {
        s1 += s2;
    }
    static T Allocate(size_t[0]) {
        return T(0);
    }
    static void ArrayFile(ofstream& file, T array, size_t[0]) {
        file << array;
    }
    template<typename V> static T Cast(V& space, size_t[0]) {
        return static_cast<T>(space);
    }
    static void Deallocate(T, size_t[0]) {
    }
    template<typename V> static void Mapping(typename Pointer<V,0>::Type& space, T[0], T*[0], size_t[0]) {
        space++;
    }
    static T Max(T space, size_t[0]) {
        return space;
    }
    static T Normalize(T& space, size_t[0], T val) {
        space /= val;
    }
    template<typename V> static T Normalize(V space, size_t[0], V val) {
        return static_cast<T>(space) / static_cast<T>(val);
    }
};
template<typename T> struct Space<T*,0> {
    static void Add(T* s1, T* s2, size_t[0]) {
        s1 += s2;
    }
    static T* Allocate(size_t[0]) {
        return new T;
    }
    static void ArrayFile(ofstream& file, T* array, size_t[0]) {
        file << *array;
    }
    template<typename V> static T* Cast(V* space, size_t[0]) {
        return static_cast<T*>(space);
    }
    static void Deallocate(T* array, size_t[0]) {
        delete array;
    }
    template<typename V> static void Mapping(typename Pointer<V*,0>::Type space, T*[0], T**[0], size_t[0]) {
        space++;
    }
    static T* Max(T* space, size_t[0]) {
        return space;
    }
    template<typename V> static T* Normalize(V space, size_t[0], V val) {
        return static_cast<T*>(space) / static_cast<T*>(val);
    }
};

#endif // SPACE_H
