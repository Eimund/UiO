/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  01.02.2015
 *
 *  c++11 compiler
 */

#ifndef ARRAY_H
#define ARRAY_H

#include "delegate.h"
#include "ref.h"
#include "pointer.h"
#include "stream.h"

using namespace std;

#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)
#define ARRAYLIST(TYPE,...) (TYPE*)(const TYPE[]){__VA_ARGS__}

template<typename C, typename T> class ArrayLength;
template<typename C, typename T> class Array;
template<typename C, typename T> class MathArray;

template<typename C, typename T> class ArrayLength {
    private: T** array;
    protected: size_t length;
    private: Delegate<C, void, T*, size_t, size_t>* init;
    public: ArrayLength(T* &array) : init(nullptr) {
        this->length = 0;
        this->array = &array;
        *this->array = new T[0];
    }
    public: ArrayLength(T* &array, size_t length) : init(nullptr) {
        this->length = 0;
        this->array = &array;
        *this->array = new T[0];
        *this = length;
    }
    public: ArrayLength(T* &array, Delegate<C, void, T*, size_t, size_t> init) : init(init) {
        this->length = 0;
        this->array = &array;
        *this->array = new T[0];
    }
    public: ArrayLength(T* &array, size_t length, Delegate<C, void, T*, size_t, size_t> init) : init(init) {
        this->length = 0;
        this->array = &array;
        *this->array = new T[0];
        *this = length;
    }
    public: size_t& operator=(const size_t& length) {
        if(length && length != this->length) {
            T* array = new T[length];
            size_t n = length < this->length ? length : this->length;
            for(size_t i = 0; i < n; i++)
                array[i] = (*this->array)[i];
            delete [] *this->array;
            *this->array = array;
            if(init != nullptr)
                (*init)(*this->array, n, length);
            return this->length = length;
        } else
            return this->length;
    }
    public: inline operator size_t () const {
        return length;
    }
    public: inline T& operator[](int i) {
        return (*array)[i];
    }
    public: inline T operator[](int i) const {
        return (*array)[i];
    }
    public: inline T* operator->() const {
        return *array;
    }
    public: inline friend Stream& operator<<(Stream& stream, const ArrayLength<C,T>& data) {
        for(size_t i = 0; i < data.length; i++)
            stream << data[i] << Stream::endl;
        return stream;
    }
    public: template<typename K, typename U> void Get(Delegate<K,U> func) {
        for(size_t i = 0; i < length; i++)
            static_cast<U&>((*array)[i]) = func();
    }
};

template<typename C, typename T> class Array : public ArrayLength<C,T> {
    friend class MathArray<C,T>;
    protected: T* array;
    public: Array() : ArrayLength<C,T>(array) {
    }
    public: Array(size_t length) : ArrayLength<C,T>(array, length) {
    }
    public: Array(Delegate<C, void, T*, size_t, size_t> init) : ArrayLength<C,T>(array, 0, init) {
    }
    public: Array(size_t length, Delegate<C, void, T*, size_t, size_t> init) : ArrayLength<C,T>(array, length, init) {
    }
    public: ~Array() {
        if(array != nullptr)
            delete [] array;
    }
    public: Array<C,T>& operator=(const Array<C,T>& other) {
        static_cast<ArrayLength<C,T>&>(*this) = static_cast<size_t>(other);
        for(size_t i = 0; i < other; i++)
            array[i] = other.array[i];
        return *this;
    }
    public: inline friend Stream& operator<<(Stream& stream, const Array<C,T>& data) {
        stream << (ArrayLength<C,T>)data;
        return stream;
    }
};

template<typename C, typename T> class MathArray : public Array<C,T> {
    public: MathArray() = default;
    public: MathArray(size_t length) : Array<C,T>(length) {
    }
    public: MathArray<C,T>& operator=(const Array<C,T>& other) {
        static_cast<ArrayLength<C,T>&>(*this) = static_cast<size_t>(other);
        for(size_t i = 0; i < other; i++)
            this->array[i] = other.array[i];
        return *this;
    }
    public: MathArray<C,T>& operator-=(T val) {
        for(size_t i = 0; i < this->length; i++)
            this->array[i] -= val;
        return *this;
    }
    public: T Sum() const {
        T sum = 0;
        for(size_t i = 0; i < this->length; i++)
            sum += this->array[i];
        return sum;
    }
};

#endif // ARRAY_H
