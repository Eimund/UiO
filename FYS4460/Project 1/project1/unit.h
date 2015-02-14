/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  04.02.2015
 *
 *  c++11 compiler
 */

#ifndef UNIT_H
#define UNIT_H

#include <cmath>
#include <string>
#include "stream.h"

using namespace std;

//----------------------------------------------------------------------------------------

template<typename V, typename... T> class Unit {
    protected: V value;
    protected: inline Unit(V value) : value(value) {
    }
    public: inline V& operator*() {
        return value;
    }
};

template<typename V, typename T> class Unit<V,T> : public Unit<V> {
    private: T* owner;
    protected: inline Unit(T* owner, V value) : Unit<V>(value), owner(owner) {
    }
    public: inline V operator*() const {
        return this->value;
    }
    public: inline T& operator=(V value) {
        this->value = value;
        return *owner;
    }
    public: inline bool operator<(const V& value) const {
        return this->value < value ? true : false;
    }
    public: inline bool operator<(const T& other) const {
        return this->value < other.value ? true : false;
    }
    public: inline bool operator>(const V& value) const {
        return this->value > value ? true : false;
    }
    public: inline bool operator>(const T& other) const {
        return this->value > other.value ? true : false;
    }
    public: inline friend Stream& operator<<(Stream& stream, const T& data) {
        stream << data.value;
        return stream;
    }
    public: inline T operator-() {
        this->value = -this->value;
        return *owner;
    }
    public: inline friend T operator+(T& v1, T& v2) {
        return *v1 + *v2;
    }
    public: inline friend T operator-(T& v1, T& v2) {
        return *v1 - *v2;
    }
    public: inline friend T operator*(T& v1, const V& v2) {
        return *v1*v2;
    }
    public: inline friend T operator/(T& v1, const V& v2) {
        return *v1/v2;
    }
    public: inline T& operator+=(const T& value) {
        this->value += *value;
        return *owner;
    }
    public: inline T& operator-=(const T& value) {
        this->value -= *value;
        return *owner;
    }
    public: inline T& operator/=(const V& value) {
        this->value /= value;
        return *owner;
    }
    public: static inline string to_string(const V& data) {
        return std::to_string(data) + " " + T::unit;
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct kg;
template<typename V> struct u;

template<typename V> struct K;

template<typename V> struct J;

template<typename V> struct m;
template<typename V> struct angstrom;

template<typename V> struct s;
template<typename V> struct ps;

template<typename V> struct J_div_K;

template<typename V> struct m_div_s;
template<typename V> struct angstrom_div_ps;

template<typename V> struct m_div_s2;

template<typename V> struct m2_div_s2;

//----------------------------------------------------------------------------------------

template<typename V> struct kg : public Unit<V,kg<V>> {
    inline kg(V value = 0) : Unit<V,kg<V>>(this, value) {
    }
    inline operator u<V>() const {
        constexpr V o = 1/1.660538921e-27;
        return this->value*o;
    }
};

template<typename V> struct u : public Unit<V,u<V>> {
    inline u(V value = 0) : Unit<V,u<V>>(this, value) {
    }
    inline operator kg<V>() const {
        constexpr V o = 1.660538921e-27;
        return this->value*o;
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct K : public Unit<V,K<V>> {
    inline K(V value = 0) : Unit<V,K<V>>(this, value) {
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct J : public Unit<V,J<V>> {
    inline J(V value = 0) : Unit<V,J<V>>(this, value) {
    }
    public: template<typename U> inline friend m2_div_s2<V> operator/(const J<V>& v1, const U& v2) {
        return (*v1)/(*v2);
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct m : public Unit<V,m<V>> {
    inline m(V value = 0) : Unit<V,m<V>>(this, value) {
    }
    inline operator angstrom<V>() const {
        constexpr V o = pow(10,10);
        return this->value*o;
    }
};

template<typename V> struct angstrom : public Unit<V,angstrom<V>> {
    static const string unit;
    inline angstrom(V value = 0) : Unit<V,angstrom<V>>(this, value) {
    }
    inline operator m<V>() const {
        constexpr V o = pow(10,-10);
        return this->value*o;
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct s : public Unit<V,s<V>> {
    static const string unit;
    inline s(V value = 0) : Unit<V,s<V>>(this, value) {
    }
    inline operator ps<V>() const {
        constexpr V o = pow(10,12);
        return this->value*o;
    }
};

template<typename V> struct ps : public Unit<V,ps<V>> {
    static const string unit;
    inline ps(V value = 0) : Unit<V,ps<V>>(this, value) {
    }
    inline operator s<V>() const {
        constexpr V o = pow(10,-12);
        return this->value*o;
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct J_div_K : public Unit<V, J_div_K<V>> {
    inline J_div_K(V value = 0) : Unit<V, J_div_K<V>>(this, value) {
    }
    public: template<typename U> inline friend J<V> operator*(const J_div_K<V>& v1, const U& v2) {
        return (*v1)*(*v2);
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct m_div_s : public Unit<V,m_div_s<V>> {
    static const string unit;
    inline m_div_s(V value = 0) : Unit<V,m_div_s<V>>(this, value) {
    }
    public: template<typename T> inline friend m<V> operator*(m_div_s<V>& v1, const T& v2) {
        return (*v1)*(*static_cast<s<V>>(v2));
    }
    inline operator angstrom_div_ps<V>() const {
        constexpr V o = pow(10,-2);
        return this->value*o;
    }
};

template<typename V> struct angstrom_div_ps : public Unit<V,angstrom_div_ps<V>> {
    static const string unit;
    inline angstrom_div_ps(V value = 0) : Unit<V,angstrom_div_ps<V>>(this, value) {
    }
    public: template<typename T> inline friend angstrom<V> operator*(angstrom_div_ps& v1, const T& v2) {
        return (*v1)*(*static_cast<ps<V>>(v2));
    }
    inline operator m_div_s<V>() const {
        constexpr V o = pow(10,2);
        return this->value*o;
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct m_div_s2 : public Unit<V,m_div_s2<V>> {
    inline m_div_s2(V value = 0) : Unit<V,m_div_s2<V>>(this, value) {
    }
    public: template<typename T> inline friend m_div_s<V> operator*(m_div_s2<V>& v1, const T& v2) {
        return (*v1)*(*static_cast<s<V>>(v2));
    }
};

//----------------------------------------------------------------------------------------

template<typename V> struct m2_div_s2 : public Unit<V,m2_div_s2<V>> {
    inline m2_div_s2(V value = 0) : Unit<V,m2_div_s2<V>>(this, value) {
    }
};

//----------------------------------------------------------------------------------------

template<typename V> const string angstrom<V>::unit = "[Å]";

template<typename V> const string s<V>::unit = "[s]";
template<typename V> const string ps<V>::unit = "[ps]";

template<typename V> const string m_div_s<V>::unit = "[m/s]";
template<typename V> const string angstrom_div_ps<V>::unit = "[Å/ps]";

//----------------------------------------------------------------------------------------

template<typename T> inline T fmod(const T& num, const T& dem) {
    return fmod(*num, *dem);
}

template<typename T> inline m_div_s<T> sqrt(const m2_div_s2<T>& value) {
    return sqrt(*value);
}

template<typename V, typename T> inline string to_string(const Unit<V,T>& val) {
    return Unit<V,T>::to_string(*val);
}

#endif // UNIT_H
