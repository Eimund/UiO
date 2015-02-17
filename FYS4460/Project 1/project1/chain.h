/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  17.02.2015
 *
 *  c++11 compiler
 */
#ifndef CHAIN_H
#define CHAIN_H

#include <cstddef>

template<class T> class Chain {
    private: T e;
    private: size_t N;
    private: Chain<T>* owner;
    private: Chain<T>* prev;
    private: Chain<T>* next;
    public: Chain() {
        N = 0;
        owner = this;
        prev = this;
        next = this;
    }
    public: Chain(const T& val) : e(val) {
        N = 0;
        owner = this;
        prev = this;
        next = this;
    }
    public: ~Chain() {
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
    public: inline operator size_t() const {
        return N;
    }
    public: inline T& operator[](const size_t& i) {
        return *static_cast<T*>(this+sizeof(*this)*i);
    }
    public: inline T operator[](const size_t& i) const {
        return *static_cast<T*>(this+sizeof(*this)*i);
    }
    public: inline void operator+=(T element) {  // Add next
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
    public: void Remove() {             // Remove next
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

#endif // CHAIN_H
