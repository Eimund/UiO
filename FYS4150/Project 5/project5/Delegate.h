/*
 *  FYS4150 - Computational Physics - Project 5
 *
 *  Written by:         Eimund Smestad
 *  Candiate number:    68
 *
 *  29.11.2014
 *
 *  c11 compiler
 */

#ifndef DELEGATE_H
#define DELEGATE_H

template<typename C, typename R, typename... P> class Delegate {
    public: typedef R (*Type)(P...);
    private: C* owner;
    private: R (C::*func)(P...);
    public: Delegate() = default;
    public: Delegate(C* owner, R (C::*func)(P...)) {
        this->owner = owner;
        this->func = func;
    }
    public: inline R operator()(P... params) {
        return (owner->*func)(params...);
    }
};
template<typename C, typename R> class Delegate<C,R,void> {
    public: typedef R (*Type)();
    private: C* owner;
    private: R (C::*func)();
    public: Delegate() {
    }
    public: Delegate(C* owner, R (C::*func)()) {
        this->owner = owner;
        this->func = func;
    }
    public: inline R operator()() {
        return (*func)();
    }
};
template<typename R, typename... P> class Delegate<void,R,P...> {
    public: typedef R (*Type)(P...);
    private: R (*func)(P...);
    public: Delegate() = default;
    public: Delegate(R (*func)(P...)) {
        this->func = func;
    }
    public: inline R operator()(P... params) {
        return (*func)(params...);
    }
};
template<typename R> class Delegate<void,R,void> {
    public: typedef R (*Type)();
    private: R (*func)();
    public: Delegate() = default;
    public: Delegate(R (*func)()) {
        this->func = func;
    }
    public: inline R operator()() {
        return (*func)();
    }
};

#endif // DELEGATE_H
