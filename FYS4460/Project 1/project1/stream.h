/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  08.02.2015
 *
 *  c++11 compiler
 */

#ifndef STREAM_H
#define STREAM_H

#include <iostream>
#include <string>

using namespace std;

class Stream {
    private: typedef Stream& (*StreamManipulator)(Stream&);
    private: typedef ofstream& (*StandardEndLine)(ofstream&);
    private: bool skip;
    private: string del;
    protected: ofstream stream;
    public: inline Stream(string del) : skip(true), del(del) {
    }
    public: inline Stream(const Stream& other) {
        this->skip = other.skip;
        this->del = other.del;
    }
    public: template<typename T> inline Stream& operator <<(const T& data) {
        if(skip) {
            stream << data;
            skip = false;
            return *this;
        }
        stream << del << data;
        return *this;
    }
    public: Stream& operator<<(StreamManipulator manip) {
        return manip(*this);
    }
    public: Stream& operator<<(StandardEndLine manip) {
        manip(stream);
        return *this;
    }
    public: inline static Stream& endl(Stream& stream) {
        stream.stream << std::endl;
        stream.skip = true;
        return stream;
    }
    public: inline void Close() {
        stream.close();
    }
    public: inline void Open(string str) {
        stream.open(str);
    }
};

#endif // STREAM_H
