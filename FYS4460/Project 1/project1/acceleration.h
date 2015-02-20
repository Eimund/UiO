/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  16.02.2015
 *
 *  c++11 compiler
 */
#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "delegate.h"

template<typename B, typename X, typename V, typename A, typename D> class Acceleration {
};
template<typename B, typename X, typename V, typename A, typename C> class Acceleration<B,X,V,A,Delegate<C,void,B&>> {
    private: Delegate<C,void,B&> f;
    public: inline void Run(B& b, X&, V&, A&) {
        f(b);
    }
};
template<typename B, typename X, typename V, typename A, typename C> struct Acceleration<B,X,V,A,Delegate<C,void,B&,A&>> {
    private: Delegate<C,void,B&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,B&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b, X&, V&, A& a) {
        f(b,a);
    }
};
template<typename B, typename X, typename V, typename A, typename C> struct Acceleration<B,X,V,A,Delegate<C,void,B&,X&,A&>> {
    private: Delegate<C,void,B&,X&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,B&,X&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b,X& x, V&, A& a) {
        f(b,x,a);
    }
};
template<typename B, typename X,typename V, typename A, typename C> struct Acceleration<B,X,V,A,Delegate<C,void,B&,V&,A&>> {
    private: Delegate<C,void,B&,V&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,B&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b,X&, V& v, A& a) {
        f(b,v,a);
    }
};
template<typename B, typename X,typename V, typename A, typename C> struct Acceleration<B,X,V,A,Delegate<C,void,B&,X&,V&,A&>> {
    private: Delegate<C,void,B&,X&,V&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,B&,X&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b, X& x, V& v, A& a) {
        f(b,x,v,a);
    }
};
template<typename X, typename V, typename A> class Acceleration<void,X,V,A,void> {
    public: static inline void Run(X&, V&, A&) {
    }
};
template<typename X, typename V, typename A, typename C> struct Acceleration<void,X,V,A,Delegate<C,void,A&>> {
    private: Delegate<C,void,A&> f;
    protected: inline Acceleration(const Delegate<C,void,A&>& f) : f(f) {
    }
    protected: inline void Run(X&, V&, A& a) {
        f(a);
    }
};
template<typename X, typename V, typename A, typename C> struct Acceleration<void,X,V,A,Delegate<C,void,X&,A&>> {
    private: Delegate<C,void,X&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,X&,A&>& f) : f(f) {
    }
    protected: inline void Run(X& x, V&, A& a) {
        f(x,a);
    }
};
template<typename X,typename V, typename A, typename C> struct Acceleration<void,X,V,A,Delegate<C,void,V&,A&>> {
    private: Delegate<C,void,V&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(X&, V& v, A& a) {
        f(v,a);
    }
};
template<typename X,typename V, typename A, typename C> struct Acceleration<void,X,V,A,Delegate<C,void,X&,V&,A&>> {
    private: Delegate<C,void,X&,V&,A&> f;
    protected: inline Acceleration(const Delegate<C,void,X&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(X& x, V& v, A& a) {
        f(x,v,a);
    }
};

#endif // ACCELERATION_H
