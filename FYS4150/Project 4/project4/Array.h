/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "Delegate.h"
#include "Type.h"

using namespace std;

#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)
#define ARRAYLIST(TYPE,...) (TYPE*)(const TYPE[]){__VA_ARGS__}

template<typename C, typename T> class ArrayLength {
    private: T** array;
    private: unsigned int length;
    private: Delegate<C, void, T*, unsigned int, unsigned int> init;
    public: ArrayLength() = default;
    public: ArrayLength(T* &array, unsigned int length, Delegate<C, void, T*, unsigned int, unsigned int> init) : init(init) {
        this->length = 0;
        this->array = &array;
        *this->array = new T[0];
        *this = length;
    }
    public: unsigned int& operator = (const unsigned int& length) {
        if(length && length != this->length) {
            T* array = new T[length];
            unsigned int n = length < this->length ? length : this->length;
            for(unsigned int i = 0; i < n; i++)
                array[i] = (*this->array)[i];
            delete [] *this->array;
            *this->array = array;
            init(*this->array, n, length);
            return this->length = length;
        } else
            return this->length;
    }
    public: inline operator unsigned int () const {
        return length;
    }
    public: inline T operator[](int i) {
        return (*array)[i];
    }
};

template<typename C, typename T> void ArrayCout(ArrayLength<C,T> array) {
    for(unsigned int i = 0; i < array; i++)
        cout << array[i] << '\t';
    cout << endl;
}
template<typename T> void ArrayCout(T* array, unsigned int n) {
    for(unsigned int i = 0; i < n; i++)
        cout << array[i] << '\t';
    cout << endl;
}
template<typename C, typename T> void ArrayToFile(ofstream& file, ArrayLength<C,T> array) {
    for(unsigned int i = 0; i < array-1; i++)
        file << array[i] << '\t';
    file << array[array-1];
}
template<typename T> void ArrayToFile(ofstream& file, T* array, unsigned int n) {
    for(unsigned int i = 0; i < n-1; i++)
        file << array[i] << '\t';
    file << array[n-1];
}

template<class T> T ArrayRelativeError(T* u, T* v, unsigned int n, T cutoff) {
    T e, error = 0;
    for(unsigned int i = 0; i < n; i++) {
        if(v[i] > cutoff || u[i] > cutoff) {
            e = abs((v[i]-u[i])/u[i]);
            if(e > error)
                error = e;
        }
    }
    return log10(error);
}

#endif // ARRAY_H
