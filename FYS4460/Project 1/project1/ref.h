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
    private: bool del;
    private: T* val;
    public: inline Ref() :del(false), val(nullptr) {
    }
    public: template<typename U> inline Ref(U val) : del(true), val(new T(val)) {
    }
    public: inline Ref(const Ref<T>& other) {
        del = false;
        val = other.val;
    }
    public: ~Ref() {
        if(del)
            delete val;
    }
    public: inline operator T& () {
        return *val;
    }
    public: inline operator const T& () const {
        return *val;
    }
    public: inline Ref<T>& operator=(T& val) {
        if(del) {
            delete this->val;
            del = false;
        }
        this->val = &val;
        return *this;
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
