/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  04.02.2015
 *
 *  c++11 compiler
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <chrono>
#include <random>
#include "unit.h"

using namespace std;

template<typename T> class Random {
    private: class _rand_ {
        private: uniform_real_distribution<T> pdf;
        private: default_random_engine rng;
        public: inline _rand_() : _rand_(0,1) {
        }
        public: inline _rand_(T min, T max) : pdf(min,max) {
            rng.seed(chrono::high_resolution_clock::now().time_since_epoch().count());
        }
        public: inline operator T() {
            return pdf(rng);
        }
    };
    public: _rand_ rand;
    public: inline Random() : rand() {
    }
    public: inline Random(T min, T max) : rand(min, max) {
    }
};

template<typename T, template<typename> class C> class Random<C<T>> {
    private: class _rand_ {
        private: Random<T> rand;
        public: inline _rand_() : Random<T>() {
        }
        public: inline operator C<T>() const {
            return static_cast<T>(rand);
        }
        public: inline T& operator*() {
            return *rand;
        }
    };
    public: _rand_ rand;
};

#endif // RANDOM_H
