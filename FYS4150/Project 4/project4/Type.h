/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */
#ifndef TYPE_H
#define TYPE_H

template<typename T> struct Type {
    static const T null;
    typedef T type;
};
template<typename T> const T Type<T>::null = 0;
template<typename T> struct Type<T*> {
    static const T* null;
    typedef T type;
};
template<typename T> const T* Type<T*>::null = nullptr;

template<typename T, typename... L> class TypeCast {
    public: typedef T type;
    protected: T value;
    protected: TypeCast(T value) : value(value) {
    }
    public: inline operator T () const {
        return value;
    }
    protected: template<typename O, typename P> struct Set {
    };
};
template<typename T, typename C> class TypeCast<T,C> : public TypeCast<T> {
    protected: TypeCast(T value) : TypeCast<T>(value) {
    }
    public: inline operator C () const {
        return this->value;
    }
    protected: template<typename O, typename P> struct Set : TypeCast<T>::template Set<O,P> {
    };
    protected: template<typename O> struct Set<O,C> {
        inline static C Value(O This, C& value) {
            This->value = value;
            return This->value;
        }
    };
};
template<typename T, typename C, typename... L>  class TypeCast<T,C,L...> : public TypeCast<T,L...> {
    protected: TypeCast(T value) : TypeCast<T,L...>(value) {
    }
    public: inline operator C () const {
        return this->value;
    }
    protected: template<typename O, typename P> struct Set : TypeCast<T,L...>:: template Set<O,P> {
    };
    protected: template<typename O> struct Set<O,C> {
        inline static C Value(O This, C& value) {
            This->value = value;
            return This->value;
        }
    };
};
#endif // TYPE_H
