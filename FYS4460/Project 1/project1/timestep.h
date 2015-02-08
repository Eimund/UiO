/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  08.02.2015
 *
 *  c++11 compiler
 */

#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <string>

using namespace std;

template<typename T, typename N> class TimeStep {
    private: T time{0};
    private: N* data;
    public: TimeStep(N& data) : data(&data) {
    }
    public: inline T& operator=(const T& time) {
        return this->time = time;
    }
    public: inline operator T&(){
        return time;
    }
    public: inline operator N*(){
        return data;
    }
    public: inline void Print() {
        *data << "Time = " + to_string(time);
    }
};

#endif // TIMESTEP_H
