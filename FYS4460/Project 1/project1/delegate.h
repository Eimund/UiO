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

template<typename... D> struct DelegateElement;
template<typename C, typename R, typename... P, typename... D> struct DelegateElement<Delegate<C,R,P...>,D...> : public DelegateElement<D...> {
    const Delegate<C,R,P...>& delegate;
    inline DelegateElement(const Delegate<C,R,P...>& delegate, const D&... delegates) : DelegateElement<D...>(delegates...), delegate(delegate) {
    }
    inline void operator()(P... params) {
        delegate(params...);
        DelegateElement<D...>::operator()(params...);
    }
};
template<typename C, typename R, typename... P> struct DelegateElement<Delegate<C,R,P...>> {
    const Delegate<C,R,P...>& delegate;
    inline DelegateElement(const Delegate<C,R,P...>& delegate) : delegate(delegate) {
    }
    inline void operator()(P... params) {
        delegate(params...);
    }
};

template<typename... P> struct DelegateArray;
template<typename... C, typename... R, typename... P> struct DelegateArray<Delegate<C,R,P...>...> : public DelegateElement<Delegate<C,R,P...>...> {
    inline DelegateArray(const Delegate<C,R,P...>&... delegates) : DelegateElement<Delegate<C,R,P...>...>(delegates...) {
    }
    inline void operator()(P... params) {
        DelegateElement<Delegate<C,R,P...>...>::operator()(params...);
    }
};

template<typename R, typename... P> inline Delegate<void,R,P...> delegate(R (*func)(P...)) {
    return Delegate<void,R,P...>(func);
}
template<typename C, typename R, typename... P> inline Delegate<C,R,P...> delegate(C& owner, R (C::*func)(P...)) {
    return Delegate<C,R,P...>(owner, func);
}
template<typename... C, typename... R, typename... P> inline DelegateArray<Delegate<C,R,P...>...> delegate_array(const Delegate<C,R,P...>&... delegates) {
    return DelegateArray<Delegate<C,R,P...>...>(delegates...);
}

#endif // DELEGATE_
