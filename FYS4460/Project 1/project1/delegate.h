/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  28.01.2015
 *
 *  c++11 compiler
 */

#ifndef DELEGATE_H
#define DELEGATE_H

template<typename C, typename R, typename... P> class Delegate {
    public: typedef R (*Type)(P...);
    private: C* owner;
    private: R (C::*func)(P...);
    public: Delegate(C& owner, R (C::*func)(P...)) : owner(&owner), func(func) {
    }
    public: inline R operator()(P... params) const {
        return (owner->*func)(params...);
    }
};
template<typename C, typename R> class Delegate<C,R,void> {
    public: typedef R (*Type)();
    private: C* owner;
    private: R (C::*func)();
    public: Delegate(C& owner, R (C::*func)()) : owner(&owner), func(func) {
    }
    public: inline R operator()() const {
        return (owner->*func)();
    }
};
template<typename R, typename... P> class Delegate<void,R,P...> {
    public: typedef R (*Type)(P...);
    private: R (*func)(P...);
    public: Delegate(R (*func)(P...)) : func(func) {
    }
    public: inline R operator()(P... params) const {
        return (*func)(params...);
    }
};
template<typename R> class Delegate<void,R,void> {
    public: typedef R (*Type)();
    private: R (*func)();
    public: Delegate(R (*func)()) : func(func) {
    }
    public: inline R operator()() const {
        return (*func)();
    }
};

#endif // DELEGATE_
