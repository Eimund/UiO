/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  20.02.2015
 *
 *  c++11 compiler
 */

#ifndef DELEGATE_H
#define DELEGATE_H

template<typename R, typename... P> class Delegate {
    private: void* owner;
    private: R (Delegate<R,P...>::*func1)(P...);
    private: R (*func2)(P...);
    public: inline Delegate() : owner(nullptr), func1(nullptr), func2(&Null<R,P...>::Function) {
    }
    public: inline Delegate(R (*func)(P...)) : owner(nullptr), func1(nullptr), func2(func) {
    }
    public: template<typename C> inline Delegate(C& owner, R (C::*func)(P...)) :
        owner(&owner), func1(reinterpret_cast<R (Delegate<R,P...>::*)(P...)>(func)), func2(nullptr) {
    }
    public: inline R operator()(P... params) const {
        if(owner)
            return (static_cast<Delegate<R,P...>*>(owner)->*func1)(params...);
        return func2(params...);
    }
    private: template<typename S, typename... Q> struct Null {
        static inline S Function(Q...) {
            return{0};
        }
    };
    private: template<typename... Q> struct Null<void,Q...> {
        static inline void Function(Q...) {
        }
    };
};

template<typename... D> struct DelegateElement;
template<typename R, typename... P, typename... D> struct DelegateElement<Delegate<R,P...>,D...> : public DelegateElement<D...> {
    const Delegate<R,P...>& delegate;
    inline DelegateElement(const Delegate<R,P...>& delegate, const D&... delegates) : DelegateElement<D...>(delegates...), delegate(delegate) {
    }
    inline void operator()(P... params) {
        delegate(params...);
        DelegateElement<D...>::operator()(params...);
    }
};
template<typename R, typename... P> struct DelegateElement<Delegate<R,P...>> {
    const Delegate<R,P...>& delegate;
    inline DelegateElement(const Delegate<R,P...>& delegate) : delegate(delegate) {
    }
    inline void operator()(P... params) {
        delegate(params...);
    }
};

template<typename... P> struct Delegating;
template<typename... R, typename... P> struct Delegating<Delegate<R,P...>...> : public DelegateElement<Delegate<R,P...>...> {
    inline Delegating(const Delegate<R,P...>&... delegates) : DelegateElement<Delegate<R,P...>...>(delegates...) {
    }
    inline void operator()(P... params) {
        DelegateElement<Delegate<R,P...>...>::operator()(params...);
    }
};

template<typename R, typename... P> inline Delegate<R,P...> delegate(R (*func)(P...)) {
    return Delegate<R,P...>(func);
}
template<typename C, typename R, typename... P> inline Delegate<R,P...> delegate(C& owner, R (C::*func)(P...)) {
    return Delegate<R,P...>(owner, func);
}
template<typename... R, typename... P> inline Delegating<Delegate<R,P...>...> delegating(const Delegate<R,P...>&... delegates) {
    return Delegating<Delegate<R,P...>...>(delegates...);
}


#endif // DELEGATE_H
