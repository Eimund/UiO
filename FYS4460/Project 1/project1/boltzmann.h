/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  04.02.2015
 *
 *  c++11 compiler
 */

#ifndef BOLTZMANN_H
#define BOLTZMANN_H

#include <cmath>
#include "unit.h"

template<typename T> class Boltzmann {
    public: K<T> temp;
    public: kg<T> mass;
    public: const J_div_K<T> kB{1.3806488e-23};
    public: Boltzmann() : temp(0), mass(1) {
    }
    public: Boltzmann(const K<T>& temp, const kg<T>& mass) : temp(temp), mass(mass) {
    }
    public: m_div_s<T> StandardDeviation() const {
        return sqrt(kB*temp/mass);
    }
};

#endif // BOLTZMANN_H
