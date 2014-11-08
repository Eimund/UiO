/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */

#ifndef PROPERTY_H
#define PROPERTY_H

#include "Delegate.h"
#include "Type.h"

enum class PropertyType {
    ReadWrite,
    ReadOnly,
    WriteOnly
};

#define PROPERTY(X,...) decltype(X)(typename decltype(X)::type(__VA_ARGS__))

template<PropertyType A, typename C, typename T, typename... L> class Property : public TypeCast<T,L...> {
    friend class Type<C>::type;
    public: typedef T type;
    private: Property() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: Property(T value) : TypeCast<T,L...>(value) {
    }
    private: inline Property& operator=(const Property& other) {
        return *this;
    }
    public: inline T operator&() {
        return this->value;
    }
    public: inline void operator=(T value) {
        this->value = value;
    }
    public: template<typename P> inline P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};
template<typename C, typename T, typename... L> class Property<PropertyType::ReadOnly, C, T, L...> : public TypeCast<T,L...> {
    friend class Type<C>::type;
    public: typedef T type;
    private: Property() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: Property(T value) : TypeCast<T,L...>(value) {
    }
    public: inline T operator&() {
        return this->value;
    }
    private: inline Property& operator=(const Property& other) {
        return *this;
    }
    private: inline void operator=(T value) {
        this->value = value;
    }
    public: template<typename P> inline P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};
template<typename C, typename T, typename... L> class Property<PropertyType::WriteOnly, C, T, L...> : protected TypeCast<T,L...> {
    friend class Type<C>::type;
    public: typedef T type;
    private: Property() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: Property(T value) : TypeCast<T,L...>(value) {
    }
    private: inline Property& operator=(const Property& other) {
        return *this;
    }
    private: inline T operator&() {
        return this->value;
    }
    public: inline void operator=(T value) {
        this->value = value;
    }
    public: template<typename P> inline P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};

template<typename C, typename T> class PropertySet {
    private: T value;
    private: Delegate<C,T,T> set;
    public: PropertySet(T value, Delegate<C,T,T> set) : value(value), set(set) {
    }
    public: T& operator = (const T& value) {
        return this->value = set(value);
    }
    public: inline operator T () const {
        return this->value;
    }
};

#endif // PROPERTY_H
