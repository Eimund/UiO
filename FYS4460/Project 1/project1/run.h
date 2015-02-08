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

#include "array.h"
#include "timestep.h"

template<typename T, typename N, typename C> void Run(TimeStep<T,N> data, Array<C,Delegate<void,void,const TimeStep<T,N>&>> func, T t_end) {
    while(data < t_end) {
        for(size_t i = 0; i < func; i++)
            func[i](data);
    }
}

#endif // RUN_H
