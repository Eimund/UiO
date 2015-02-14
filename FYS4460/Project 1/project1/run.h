/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  08.02.2015
 *
 *  c++11 compiler
 */

#ifndef RUN_H
#define RUN_H

#include "delegate.h"
#include "timestep.h"

template<typename T, typename N, typename... C, typename... R> void Run(TimeStep<T,N>& data, T t_end, const Delegate<C,R,TimeStep<T,N>&>&... f) {
    auto func = delegate_array(f...);
    while(static_cast<T&>(data) < t_end) {
        data.Print();
        func(data);
    }
    data.Save();
}

#endif // RUN_H
