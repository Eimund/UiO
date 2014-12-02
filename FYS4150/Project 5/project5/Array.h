/*
 *  FYS4150 - Computational Physics - Project 5
 *
 *  Written by:         Eimund Smestad
 *  Candiate number:    68
 *
 *  29.11.2014
 *
 *  c++11 compiler
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

template<class T> struct Chain {
    unsigned int N;
    T e;
    Chain<T>* owner;
    Chain<T>* prev;
    Chain<T>* next;
    Chain() {
        N = 0;
        owner = this;
        prev = this;
        next = this;
    }
    ~Chain() {
        if(N) {
            Chain<T>* p = this;
            for(int i = 1; i < N; i++)
                p = p->next;
            for(int i = 0; i < N; i++) {
                delete p->next;
                p = p->prev;
            }
        }
    }
    void Add(T element) {       // Add next
        if(next == this) {
            next = new Chain<T>;
            next->next = next;
        }
        else {
            next->prev = new Chain<T>;
            next->prev->next = next;
            next = next->prev;
        }
        owner->N++;
        next->owner = owner;
        next->prev = this;
        next->e = element;
    }
    void Remove() {             // Remove next
        if(next->next == next) {
            delete next;
            next = this;
        } else {
            next = next->next;
            delete next->prev;
            next->prev = this;
        }
        owner->N--;
    }
};

template<typename T, unsigned int D> struct Pointer : Pointer<T*,D-1> {
};
template<typename T> struct Pointer<T,0> {
    typedef T Type;
};

template<typename T, unsigned int D> struct Vector {
    T e[D];
    operator T*() {
        return e;
    }
    Vector & operator= (const T other[D]) {
        for(int i = 0; i < D; i++)
            e[i] = other[i];
    }
};

template<typename T, unsigned int D> struct Space {
    static void Add(typename Pointer<T,D>::Type s1, typename Pointer<T,D>::Type s2, unsigned int n[D]) {
        for(int i = 0; i < n[0]; i++)
            Space<T,D-1>::Add(s1[i], s2[i], &n[1]);
    }
    static typename Pointer<T,D>::Type Allocate(unsigned int n[D]) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(int i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::Allocate(&n[1]);
        return array;
    }
    static void ArrayFile(ofstream& file, typename Pointer<T,D>::Type array, unsigned int n[D]) {
        for(int i = 0; i < n[0]; i++) {
            if(D == 1 && i)
                file << '\t';
            else if(D > 1 && i)
                file << endl;
            Space<T,D-1>::ArrayFile(file, array[i], &n[1]);
        }
    }
    template<typename V> static typename Pointer<T,D>::Type Cast(typename Pointer<V,D>::Type space, unsigned int n[D]) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(int i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::Cast(space[i], &n[1]);
        return array;
    }
    static void Deallocate(typename Pointer<T,D>::Type array, unsigned int n[D]) {
        for(int i = 0; i < n[0]; i++)
            Space<T,D-1>::Deallocate(array[i], &n[1]);
        delete [] array;
    }
    static void DeRange(T** range) {
        for(int i = 0; i < D; i++)
            delete [] range[i];
        delete [] range;
    }
    template<typename V> static typename Pointer<V,D>::Type Map(Chain<Vector<T,D>>& chain, T* range[D], unsigned int n[D]) {
        auto space = Space<V,D>::Allocate(n);
        auto p = chain.owner;
        for(unsigned int i = 0; i < chain.N; i++) {
            p = p->next;
            Space<T,D>::Mapping<V>(space, p->e.e, range, n);
        }
        return space;
    }
    template<typename V> static void Map(typename Pointer<V,D>::Type& space, Chain<Vector<T,D>>& chain, T* range[D], unsigned int n[D]) {
        auto p = chain.owner;
        for(unsigned int i = 0; i < chain.N; i++) {
            p = p->next;
            Space<T,D>::Mapping<V>(space, p->e.e, range, n);
        }
    }
    template<typename V> static void Mapping(typename Pointer<V,D>::Type space, T element[D], T* range[D], unsigned int n[D]) {
        unsigned int m = n[0]-1;
        if(element[0] >= range[0][0] || element[0] <= range[0][m]) {
            for(unsigned int i = 0; i < m; i++) {
                if(element[0] < (range[0][i]+range[0][i+1])/2) {
                    Space<T,D-1>::template Mapping<V>(space[i], &element[1], &range[1], &n[1]);
                    return;
                }
            }
            Space<T,D-1>::template Mapping<V>(space[m], &element[1], &range[1], &n[1]);
        }
    }
    static T Max(typename Pointer<T,D>::Type space, unsigned int n[D]) {
        T val, max = 0;
        for(unsigned int i = 0; i < n[0]; i++) {
            val = Space<T,D-1>::Max(space[i], &n[1]);
            if(val > max)
                max = val;
        }
        return max;
    }
    static void Normalize(typename Pointer<T,D>::Type space, unsigned int n[D], T val) {
        for(int i = 0; i < n[0]; i++)
            Space<T,D-1>::Normalize(space[i], &n[1], val);
    }
    template<typename V> static typename Pointer<T,D>::Type Normalize(typename Pointer<V,D>::Type space, unsigned int n[D], V val) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(int i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::template Normalize<V>(space[i], &n[1], val);
        return array;
    }
    static T** Range(T lower[D], T upper[D], unsigned int n[D]) {
        T** array = new T*[D];
        for(int i = 0; i < D; i++) {
            T step = (upper[i]-lower[i])/(n[i]-1);
            array[i] = new T[n[i]];
            array[i][0] = lower[i];
            for(int j = 1; j < n[i]; j++)
                array[i][j] = array[i][j-1] + step;
        }
        return array;
    }
    static void RangeFile(ofstream& file, T* range[D], unsigned int n[D]) {
        for(int i = 0; i < D; i++) {
            if(i)
                file << endl;
            Space<T,1>::ArrayFile(file, range[i], &n[i]);
        }
    }
    static void ToFile(ofstream& file, T* range[D], typename Pointer<T,D>::Type array, unsigned int n[D]) {
        file << D << endl;
        RangeFile(file, range, n);
        file << endl;
        ArrayFile(file, array, n);
    }
};
template<typename T> struct Space<T,0> {
    static void Add(T& s1, T& s2, unsigned int[0]) {
        s1 += s2;
    }
    static T Allocate(unsigned int[0]) {
        return T(0);
    }
    static void ArrayFile(ofstream& file, T array, unsigned int[0]) {
        file << array;
    }
    template<typename V> static T Cast(V& space, unsigned int[0]) {
        return static_cast<T>(space);
    }
    static void Deallocate(T, unsigned int[0]) {
    }
    template<typename V> static void Mapping(typename Pointer<V,0>::Type& space, T[0], T*[0], unsigned int[0]) {
        space++;
    }
    static T Max(T space, unsigned int[0]) {
        return space;
    }
    static T Normalize(T& space, unsigned int[0], T val) {
        space /= val;
    }
    template<typename V> static T Normalize(V space, unsigned int[0], V val) {
        return static_cast<T>(space) / static_cast<T>(val);
    }
};
template<typename T> struct Space<T*,0> {
    static void Add(T* s1, T* s2, unsigned int[0]) {
        s1 += s2;
    }
    static T* Allocate(unsigned int[0]) {
        return new T;
    }
    static void ArrayFile(ofstream& file, T* array, unsigned int[0]) {
        file << *array;
    }
    template<typename V> static T* Cast(V* space, unsigned int[0]) {
        return static_cast<T*>(space);
    }
    static void Deallocate(T* array, unsigned int[0]) {
        delete array;
    }
    template<typename V> static void Mapping(typename Pointer<V*,0>::Type space, T*[0], T**[0], unsigned int[0]) {
        space++;
    }
    static T* Max(T* space, unsigned int[0]) {
        return space;
    }
    template<typename V> static T* Normalize(V space, unsigned int[0], V val) {
        return static_cast<T*>(space) / static_cast<T*>(val);
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
template<typename T> unsigned int ArrayFind(T* x, T dx, unsigned int i, unsigned int n) {
    T xj = x[i] + dx;
    if(dx > 0) {
        for(int j = i; j < n; j++) {
            if(xj <= x[j]) {
                if(j == 0)
                    return 0;
                if(abs(xj-x[j]) < abs(xj-x[j-1]))
                    return j;
                return j-1;
            }
        }
    }
    for(int j = i; j; j--) {
        if(xj >= x[j]) {
            if(j == n-1)
                return n-1;
            if(abs(xj-x[j]) < abs(xj-x[j+1]))
                return j;
            return j+1;
        }
    }
    return 0;
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
