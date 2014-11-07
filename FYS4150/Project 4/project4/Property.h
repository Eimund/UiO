/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */

#ifndef PROPERTY_H
#define PROPERTY_H

template<typename T, typename C, typename R, typename P> class Property {
    Delegate<C,void,P> Set;
    Delegate<C,R,void> Get;

};

#endif // PROPERTY_H
