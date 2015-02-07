/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  07.02.2015
 *
 *  c++11 compiler
 */

#ifndef TEMPLATE_H
#define TEMPLATE_H

template<typename T> struct RemoveTemplate {
    typedef T type;
};
template<typename T, template<typename> class C> struct RemoveTemplate<C<T>> : RemoveTemplate<T> {

};

#endif // TEMPLATE_H
