/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  28.01.2015
 *
 *  c++11 compiler
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <memory>
#include <string>
#include "array.h"

using namespace std;

template<typename T, size_t D> struct Vector {
    T e[D];
    Vector() : Vector(static_cast<T>(0)) {
    }
    Vector(const T& val) {
        for(size_t i = 0; i < D; i++)
            e[i] = val;
    }
    Vector(const T val[D]) : e(val) {
    }
    template<typename... P> Vector(P... val) : e{val...} {
    }
    inline operator T*() {
        return e;
    }
    inline Vector& operator= (const T other[D]) {
        for(size_t i = 0; i < D; i++)
            e[i] = other[i];
    }
    inline T& operator[](size_t i) {
        return e[i];
    }
    inline T operator[](size_t i) const {
        return e[i];
    }
    inline operator size_t() const {
        return D;
    }
    inline friend Stream& operator<<(Stream& stream, const Vector<T,D>& data) {
        for(size_t i = 0; i < D; i++)
            stream << data[i];
        return stream;
    }
};

template<typename T, size_t D, size_t N> struct Vectors {
    Vector<T,D> v[N];
    Vectors(T val = 0) {
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < D; j++)
                v[i].e[j] = val;
        }
    }
    inline Vector<T,D>& operator[](size_t i) {
        return v[i];
    }
    template<typename V> inline Array<Ref<T>> Get(ArrayLength<V> other, size_t j, size_t k) {
        Array<Ref<T>> data(other);
        for(size_t i = 0; i < other; i++)
            data[i] = other[i][j][k];
        return data;
    }
};

template<typename T, size_t D, size_t N> struct VectorLabel : Vectors<T,D,N> {
    string l;
    VectorLabel(T val = 0) : Vectors<T,D,N>(val) {
    }
    inline operator string() const {
        return l;
    }
    inline VectorLabel<T,D,N>& operator=(const string l) {
        this->l = l;
        return *this;
    }
};

#endif // VECTOR_H
