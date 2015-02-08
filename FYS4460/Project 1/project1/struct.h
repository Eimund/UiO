/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  08.02.2015
 *
 *  c++11 compiler
 */

#ifndef STRUCT_H
#define STRUCT_H

#include "stream.h"

template<size_t N, typename T> struct StructData {
    T data;
    inline friend Stream& operator<<(Stream& stream, const StructData<N,T>& data) {
        stream << data.data;
        return stream;
    }
};
template<size_t N, typename T, typename... P> struct StructMore : StructMore<N+1,P...>, StructData<N,T> {
    inline friend Stream& operator<<(Stream& stream, const StructMore<N,T,P...>& data) {
        stream << static_cast<StructData<N,T>>(data);
        stream << static_cast<StructMore<N+1,P...>>(data);
        return stream;
    }
};
template<size_t N, typename T> struct StructMore<N,T> : StructData<N,T> {
    inline friend Stream& operator<<(Stream& stream, const StructMore<N,T>& data) {
        stream << static_cast<StructData<N,T>>(data);
        return stream;
    }
};
template<typename... P> struct Struct : StructMore<0,P...> {
    inline friend Stream& operator<<(Stream& stream, const Struct<P...>& data) {
        stream << static_cast<StructMore<0,P...>>(data);
        return stream;
    }
};

#endif // STRUCT_H
