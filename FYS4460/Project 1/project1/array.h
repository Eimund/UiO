/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  20.02.2015
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

template<typename T> class ArrayLength;
template<typename T> class Array;
template<typename T> class MathArray;

template<typename T> class ArrayLength {
    private: T** array;
    protected: size_t length{0};
    public: Delegate<void,T*,size_t,size_t> init;
    public: ArrayLength(T* &array) : array(&array) {
        *this->array = new T[0];
    }
    public: ArrayLength(T* &array, size_t length) : array(&array) {
        *this->array = new T[0];
        *this = length;
    }
    public: ArrayLength(T* &array, const Delegate<void,T*,size_t,size_t>& init) : array(&array), init(init) {
        *this->array = new T[0];
    }
    public: ArrayLength(T* &array, size_t length, const Delegate<void,T*,size_t,size_t>& init) : array(&array), init(init) {
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
            init(*this->array, n, length);
        }
        return this->length = length;
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
    public: inline friend Stream& operator<<(Stream& stream, const ArrayLength<T>& data) {
        for(size_t i = 0; i < data.length; i++)
            stream << data[i] << Stream::endl;
        return stream;
    }
    public: template<typename U> void Get(Delegate<U> func) {
        for(size_t i = 0; i < length; i++)
            static_cast<U&>((*array)[i]) = func();
    }
    public: void SetInit(Delegate<void,T*,size_t,size_t>& init) {
        this->init = &init;
    }
};

template<typename T> class Array : public ArrayLength<T> {
    friend class MathArray<T>;
    protected: T* array;
    public: Array() : ArrayLength<T>(array) {
    }
    public: Array(size_t length) : ArrayLength<T>(array, length) {
    }
    public: Array(Delegate<void,T*,size_t,size_t>& init) : ArrayLength<T>(array, 0, init) {
    }
    public: Array(size_t length, const Delegate<void,T*,size_t,size_t>& init) : ArrayLength<T>(array, length, init) {
    }
    public: ~Array() {
        if(array != nullptr)
            delete [] array;
    }
    public: Array(const Array<T>& other) : ArrayLength<T>(array, other.length, other.init) {
        *this = other;
    }
    public: inline size_t& operator=(const size_t& length) {
        return static_cast<ArrayLength<T>&>(*this) = length;
    }
    public: Array<T>& operator=(const Array<T>& other) {
        *this = static_cast<size_t>(other);
        this->init = other.init;
        for(size_t i = 0; i < other; i++)
            array[i] = other.array[i];
        return *this;
    }
    public: inline friend Stream& operator<<(Stream& stream, const Array<T>& data) {
        stream << (ArrayLength<T>)data;
        return stream;
    }
};

template<typename T> class MathArray : public Array<T> {
    public: MathArray() = default;
    public: MathArray(size_t length) : Array<T>(length) {
    }
    public: MathArray<T>& operator=(const Array<T>& other) {
        static_cast<ArrayLength<T>&>(*this) = static_cast<size_t>(other);
        for(size_t i = 0; i < other; i++)
            this->array[i] = other.array[i];
        return *this;
    }
    public: MathArray<T>& operator-=(const T& val) {
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
