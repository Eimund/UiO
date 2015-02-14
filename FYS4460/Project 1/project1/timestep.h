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
    private: size_t i{0};
    private: size_t n;
    private: T time{0};
    public: N* data;
    public: inline TimeStep(N& data, size_t n = 1) : n(n), data(&data) {
    }
    public: inline T& operator=(const T& time) {
        return this->time = time;
    }
    public: inline operator size_t() {
        return static_cast<size_t>(*data);
    }
    public: inline operator T&(){
        return time;
    }
    public: inline operator N*(){
        return data;
    }
    public: inline void Close() {
        data->Close();
    }
    public: inline void Print() {
        if(i == 0)
            Save();
        if(++i >= n)
            i = 0;
    }
    public: inline void Open(string str) {
        data->Open(str);
        i = 0;
    }
    public: inline void Save() {
        *data << "Time = " + to_string(time);
    }
};

#endif // TIMESTEP_H
