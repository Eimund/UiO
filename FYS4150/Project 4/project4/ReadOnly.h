/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */
#ifndef READONLY_H
#define READONLY_H

#include "Type.h"

template<typename C, typename T, typename... L> class ReadOnly : public TypeCast<T,L...> {
    friend class Type<C>::type;
    private: ReadOnly() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: ReadOnly(T value) : TypeCast<T,L...>(value) {
    }
    private: ReadOnly& operator=(const ReadOnly& other) {
        return *this;
    }
    public: template<typename P> P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};

#endif // READONLY_H
