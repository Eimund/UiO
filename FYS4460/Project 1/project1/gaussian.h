#ifndef GAUSSIAN_H
#define GAUSSIAN_H

/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  31.01.2015
 *
 *  c++11 compiler
 */

#include <cmath>
#include "random.h"
#include "template.h"
#include "vector.h"

using namespace std;

template<typename T> class Gaussian : public Random<typename RemoveTemplate<T>::type> {
    public: T mu;
    public: T sigma;
    public: Gaussian(T sigma) : Gaussian(0, sigma) {
    }
    public: Gaussian(const T& mu, const T& sigma) : Random<typename RemoveTemplate<T>::type>(-1,1), mu(mu), sigma(sigma) {
    }
    public: inline T Distribution() {
        typename RemoveTemplate<T>::type val = this->rand;
        if(val < 0)
            return mu-sigma*sqrt(-2*log(abs(val)));
        return mu+sigma*sqrt(-2*log(abs(val)));
    }
    public: inline T Distribution_mu_0() {
        typename RemoveTemplate<T>::type val = this->rand;
        if(val < 0)
            return -sigma*sqrt(-2*log(abs(val)));
        return sigma*sqrt(-2*log(abs(val)));
    }
};

#endif // GAUSSIAN_H
