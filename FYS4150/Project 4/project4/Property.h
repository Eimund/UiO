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

/*template<typename T, typename C, typename R, typename P> class Property {
    Delegate<C,void,P> Set;
    Delegate<C,R,void> Get;

};*/

#include "Type.h"

enum class PropertyType {
    ReadWrite,
    ReadOnly,
    WriteOnly
};

template<PropertyType A, typename C, typename T, typename... L> class Property : public TypeCast<T,L...> {
    friend class Type<C>::type;
    private: Property() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: Property(T value) : TypeCast<T,L...>(value) {
    }
    private: Property& operator=(const Property& other) {
        return *this;
    }
    public: T operator&() {
        return this->value;
    }
    public: void operator=(T value) {
        this->value = value;
    }
    public: template<typename P> P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};
template<typename C, typename T, typename... L> class Property<PropertyType::ReadOnly, C, T, L...> : public TypeCast<T,L...> {
    friend class Type<C>::type;
    private: Property() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: Property(T value) : TypeCast<T,L...>(value) {
    }
    public: T operator&() {
        return this->value;
    }
    private: Property& operator=(const Property& other) {
        return *this;
    }
    private: void operator=(T value) {
        this->value = value;
    }
    public: template<typename P> P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};
template<typename C, typename T, typename... L> class Property<PropertyType::WriteOnly, C, T, L...> : protected TypeCast<T,L...> {
    friend class Type<C>::type;
    private: Property() : TypeCast<T,L...>(Type<T>::null) {
    }
    private: Property(T value) : TypeCast<T,L...>(value) {
    }
    private: Property& operator=(const Property& other) {
        return *this;
    }
    private: T operator&() {
        return this->value;
    }
    public: void operator=(T value) {
        this->value = value;
    }
    public: template<typename P> P operator = (P value) {
        return TypeCast<T,L...>::template Set<decltype(this),P>::Value(this, value);
    }
};

#endif // PROPERTY_H
