#ifndef DELEGATE_H
#define DELEGATE_H

template<typename C, typename R, typename... P> class Delegate {
    public: typedef R (*Type)(P...);
    private: C* owner;
    private: R (C::*func)(P...);
    public: Delegate(C* owner, R (C::*func)(P...)) {
        this->owner = owner;
        this->func = func;
    }
    public: R operator()(P... params) {
        return (owner->*func)(params...);
    }
};

template<typename C, typename R> class Delegate<C,R,void> {
    public: typedef R (*Type)();
    private: C* owner;
    private: R (C::*func)();
    public: Delegate(C* owner, R (C::*func)()) {
        this->owner = owner;
        this->func = func;
    }
    public: R operator()() {
        return (owner->*func)();
    }
};

#endif // DELEGATE_H
