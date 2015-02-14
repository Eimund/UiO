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
#include <fstream>
#include <string>

using namespace std;

class Stream {
    private: typedef Stream& (*StreamManipulator)(Stream&);
    private: typedef ofstream& (*StandardEndLine)(ofstream&);
    protected: bool first;
    private: bool newline;
    private: string del;
    protected: ofstream stream;
    public: inline Stream(string del) : first(true), newline(false), del(del) {
    }
    public: inline Stream(const Stream& other) {
        this->first = other.first;
        this->newline = other.newline;
        this->del = other.del;
    }
    public: template<typename T> inline Stream& operator <<(const T& data) {
        if(first) {
            stream << data;
            first = false;
            newline = false;
            return *this;
        }
        if(newline) {
            stream << std::endl << data;
            newline = false;
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
        stream.newline = true;
        return stream;
    }
    public: inline void Close() {
        stream.close();
    }
    public: inline void Open(string str) {
        stream.open(str);
        first = true;
        newline = false;
    }
};

#endif // STREAM_H
