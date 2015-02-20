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
    private: size_t N{0};
    private: size_t ilast{0};
    private: Chain<T>* owner;
    private: Chain<T>* prev;
    private: Chain<T>* next;
    private: Chain<T>* last;
    public: inline Chain()
        : owner(this), prev(this), next(this), last(this) {
    }
    public: inline Chain(const T& val)
        : e(val), owner(this), prev(this), next(this), last(this) {
    }
    public: inline Chain(const T& val, Chain<T>* owner, Chain<T>* prev)
        : e(val), owner(owner), prev(prev), next(this), last(this) {
    }
    public: inline Chain(const T& val, Chain<T>* owner, Chain<T>* prev, Chain<T>* next)
        : e(val), owner(owner), prev(prev), next(next), last(this) {
    }
    public: inline Chain(const Chain<T>& other) {
        *this = other;
    }
    public: inline ~Chain() {
        Free();
    }
    public: inline Chain<T>& operator=(const Chain<T>& other) {
        if(&other == other.owner) {
            Free();
            this->owner = this;
            this->next = this;
            this->prev = this;
            this->ilast = 0;
            this->last = this;
            auto p = other.next;
            for(size_t i = 0; i < other.N; i++) {
                *this += *p;
                p = p->next;
            }
            /*this->last = this->next;
            this->ilast = 0;
            this->N = other.N;
            auto p = this->next;
            for(size_t i = 0; i < this->N; i++) {
                p->owner = this;
                p = p->next;
            }*/

        }
        this->e = other.e;
        return *this;
    }
    public: inline operator size_t() const {
        return N;
    }
    public: inline T& operator[](const size_t& i) {
        return Index(i).e;
    }
    public: inline T operator[](const size_t& i) const {
        return Index(i).e;
    }
    public: inline Chain<T>& operator+=(const T& element) {         // Add next LIFO
        if(next == this)
            next = new Chain<T>(element, owner, this);
        else {
            next->prev = new Chain<T>(element, owner, this, next);
            next = next->prev;
        }
        last = owner->next;
        ilast = 0;
        owner->N++;
        return *this;
    }
    public: inline Chain<T>& operator+=(const Chain<T>& element) {  // Add next FIFO
        auto p = N ? Index(N-1) : *this;
        p += element.e;
        return *this;
    }
    public: inline void Free() {
        if(N) {
            auto p = this;
            for(size_t i = 1; i < N; i++)
                p = p->next;
            for(size_t i = 0; i < N; i++) {
                delete p->next;
                p = p->prev;
            }
            N = 0;
        }
    }
    private: inline Chain<T>& Index(const size_t& i) {
        if(i <= ilast) {
            for(size_t j = ilast-i; j; j--)
                last = last->prev;
        } else {
            for(size_t j = i-ilast; j; j--)
                last = last->next;
        }
        ilast = i;
        return *last;
    }
    public: inline void Remove() {                      // Remove next
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
