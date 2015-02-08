/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  28.01.2015
 *
 *  c++11 compiler
 */

#ifndef POINTER_H
#define POINTER_H

template<typename T, size_t D> struct Pointer : Pointer<T*,D-1> {
};
template<typename T> struct Pointer<T,0> {
    typedef T Type;
};

#endif // POINTER_H
