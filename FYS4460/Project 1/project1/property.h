/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  15.02.2015
 *
 *  c++11 compiler
 */

#ifndef PROPERTY_H
#define PROPERTY_H

#include "delegate.h"

template<typename T> class Property {
    private: T data;
    public: Delegate<void,T&,const T&>* func;
    public: inline Property() : func(nullptr) {
    }
    public: inline Property(const T& data) : data(data), func(nullptr) {
    }
    public: inline T& operator=(const T& data) {
        if(func == nullptr)
            return this->data = data;
        (*func)(this->data, data);
        return this->data;
    }
    public: inline operator T() const {
        return data;
    }
};

#endif // PROPERTY_H
