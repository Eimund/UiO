/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  14.02.2015
 *
 *  c++11 compiler
 */
#ifndef LJPOTENTIAL_H
#define LJPOTENTIAL_H

#include "unit.h"

template<typename E, typename X, typename F> class LJ_Potential {
    private: E epsilon;
    private: X sigma;
    private: F force;
    LJ_Potential(E epsilon, X sigma) : epsilon(epsilon), sigma(sigma), force(4*epsilon/sigma) {
    }

};


#endif // LJPOTENTIAL_H
