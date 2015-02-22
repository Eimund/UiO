/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  28.01.2015
 *
 *  c++11 compiler
 */

#ifndef REF_H
#define REF_H

template<typename T> class Ref {
    private: T ref;
    private: T* val;
    public: inline Ref() : ref(0), val(&ref) {
    }
    public: inline Ref(T& val) : ref(0), val(&val) {
    }
    public: template<typename U> inline Ref(const U& val) : ref(val), val(&ref) {
    }
    public: inline Ref(const Ref<T>& other) {
        *this = other;
    }
    public: inline operator T& () {
        return *val;
    }
    public: inline operator const T& () const {
        return *val;
    }
    public: inline Ref<T>& operator=(const Ref<T>& other) {
        this->ref = other.ref;
        if(other.val == &other.ref)
            this->val = &this->ref;
        else
            this->val = other.val;
        return *this;
    }
    public: inline Ref<T>& operator=(T& val) {
        this->val = &val;
        return *this;
    }
    public: inline auto operator[](const size_t i) -> decltype((*val)[0])& {
        return (*val)[i];
    }
    public: template<typename V> inline Ref<T>& operator+=(const V& val) {
        *this->val += val;
        return *this;
    }
    public: template<typename V> inline Ref<T>& operator-=(const V& val) {
        *this->val -= val;
        return *this;
    }
};

#endif // REF_H
