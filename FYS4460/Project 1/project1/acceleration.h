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

template<typename X, typename V, typename A, typename D> class Acceleration {
};
/*template<typename B, typename X, typename V, typename A> class Acceleration<B,X,V,A,Delegate<void,B&>> {
    private: Delegate<void,B&> f;
    public: inline void Run(B& b, X&, V&, A&) {
        f(b);
    }
};
template<typename B, typename X, typename V, typename A> struct Acceleration<B,X,V,A,Delegate<void,B&,A&>> {
    private: Delegate<void,B&,A&> f;
    protected: inline Acceleration(const Delegate<void,B&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b, X&, V&, A& a) {
        f(b,a);
    }
};
template<typename B, typename X, typename V, typename A> struct Acceleration<B,X,V,A,Delegate<void,B&,X&,A&>> {
    private: Delegate<void,B&,X&,A&> f;
    protected: inline Acceleration(const Delegate<void,B&,X&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b,X& x, V&, A& a) {
        f(b,x,a);
    }
};
template<typename B, typename X,typename V, typename A> struct Acceleration<B,X,V,A,Delegate<void,B&,V&,A&>> {
    private: Delegate<void,B&,V&,A&> f;
    protected: inline Acceleration(const Delegate<void,B&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b,X&, V& v, A& a) {
        f(b,v,a);
    }
};
template<typename B, typename X,typename V, typename A> struct Acceleration<B,X,V,A,Delegate<void,B&,X&,V&,A&>> {
    private: Delegate<void,B&,X&,V&,A&> f;
    protected: inline Acceleration(const Delegate<void,B&,X&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(B& b, X& x, V& v, A& a) {
        f(b,x,v,a);
    }
};*/
template<typename X, typename V, typename A> class Acceleration<X,V,A,void> {
    public: static inline void Run(X&, V&, A&) {
    }
};
template<typename X, typename V, typename A> struct Acceleration<X,V,A,Delegate<void,A&>> {
    private: Delegate<void,A&> f;
    protected: inline Acceleration(const Delegate<void,A&>& f) : f(f) {
    }
    protected: inline void Run(X&, V&, A& a) {
        f(a);
    }
};
template<typename X, typename V, typename A> struct Acceleration<X,V,A,Delegate<void,X&,A&>> {
    private: Delegate<void,X&,A&> f;
    protected: inline Acceleration(const Delegate<void,X&,A&>& f) : f(f) {
    }
    protected: inline void Run(X& x, V&, A& a) {
        f(x,a);
    }
};
template<typename X,typename V, typename A> struct Acceleration<X,V,A,Delegate<void,V&,A&>> {
    private: Delegate<void,V&,A&> f;
    protected: inline Acceleration(const Delegate<void,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(X&, V& v, A& a) {
        f(v,a);
    }
};
template<typename X,typename V, typename A> struct Acceleration<X,V,A,Delegate<void,X&,V&,A&>> {
    private: Delegate<void,X&,V&,A&> f;
    protected: inline Acceleration(const Delegate<void,X&,V&,A&>& f) : f(f) {
    }
    protected: inline void Run(X& x, V& v, A& a) {
        f(x,v,a);
    }
};

#endif // ACCELERATION_H
