/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  08.02.2015
 *
 *  c++11 compiler
 */

#ifndef VMD_H
#define VMD_H

#include "stream.h"
#include "struct.h"
#include "unit.h"
#include "vector.h"

template<typename T, typename... P> struct VMD_Data : public Struct<Vector<T,3>,P...> {
    string label;
    public: inline friend Stream& operator<<(Stream& stream, const VMD_Data<T,P...>& data) {
        stream << data.label << static_cast<Struct<Vector<T,3>,P...>>(data);
        return stream;
    }
};

template<typename T, typename... P> class VMD : public Stream {
    protected: string comment{""};
    protected: Array<void,VMD_Data<T,P...>> data;
    public: inline VMD(string comment = "") : Stream(" "), comment(comment) {
    }
    public: inline VMD& operator<<(const string& data) {
        if(!this->first) {
            this->stream << std::endl;
            first = true;
        }
        this->stream << (size_t)this->data << std::endl;
        this->stream << data << ": " << this->comment << std:: endl;
        static_cast<Stream&>(*this) << this->data;
        return *this;
    }
    public: inline Vector<T,3>& operator[](size_t i) {
        return static_cast<StructData<0,Vector<T,3>>&>(data[i]).data;
    }
    public: inline Vector<T,3> operator[](size_t i) const {
        return static_cast<StructData<0,Vector<T,3>>&>(data[i]).data;
    }
    public: inline string& label(size_t i) {
        return data[i].label;
    }
    public: inline operator size_t() const {
        return data;
    }
    public: void LatticeCenteredCubic(string e[4], T b, size_t N[3]) {
        data = 4*N[0]*N[1]*N[2];
        for(size_t i = 0, j = 0; j < N[0]; j++) {
            for(size_t k = 0; k < N[1]; k++) {
                for(size_t l = 0; l < N[2]; l++) {
                    label(i) = e[0];             // Particle positions
                    x(i)[0] = b*l;
                    x(i)[1] = b*k;
                    x(i)[2] = b*j;
                    i++;
                    label(i) = e[1];
                    x(i)[0] = b*(l+0.5);
                    x(i)[1] = b*(k+0.5);
                    x(i)[2] = b*j;
                    i++;
                    label(i) = e[2];
                    x(i)[0] = b*l;
                    x(i)[1] = b*(k+0.5);
                    x(i)[2] = b*(j+0.5);
                    i++;
                    label(i) = e[3];
                    x(i)[0] = b*(l+0.5);
                    x(i)[1] = b*k;
                    x(i)[2] = b*(j+0.5);
                    i++;
                }
            }
        }
    }
    public: inline Vector<T,3>& x(size_t i) {
        return static_cast<StructData<0,Vector<T,3>>&>(data[i]).data;
    }
    public: inline Vector<T,3> x(size_t i) const {
        return static_cast<StructData<0,Vector<T,3>>&>(data[i]).data;
    }
};

template<typename T1, typename T2, typename... P> class VMD_Vel : public VMD<T1,Vector<T2,3>,P...> {
    public: inline VMD_Vel(string comment = "") : VMD<T1,Vector<T2,3>,P...>(comment) {
    }
    public: inline VMD_Vel& operator=(const VMD_Vel& other) {
        this->comment = other.comment;
        this->data = other.data;
        return *this;
    }
    public: inline Vector<T2,3>& operator[](size_t i) {
        return static_cast<StructData<1,Vector<T2,3>>&>(this->data[i]).data;
    }
    public: inline Vector<T2,3> operator[](size_t i) const {
        return static_cast<StructData<1,Vector<T2,3>>&>(this->data[i]).data;
    }
    public: template<typename C, typename U> void SetVelocityDistribution(Delegate<C,U>& distribution) {
        Vector<T2,3> sum(0);
        for(size_t i = 0; i < this->data; i++) {    // Get velocity distribution
            for(size_t j = 0; j < 3; j++) {
                v(i)[j] = distribution();
                sum[j] += v(i)[j];
            }
        }

        for(size_t j = 0; j < 3; j++)               // Mean value
            sum[j] /= this->data;

        for(size_t i = 0; i < this->data; i++) {    // Remove initial linear momentum
            for(size_t j = 0; j < 3; j++)
                v(i)[j] -= sum[j];
        }
    }
    public: inline Vector<T2,3>& v(size_t i) {
        return static_cast<StructData<1,Vector<T2,3>>&>(this->data[i]).data;
    }
    public: inline Vector<T2,3> v(size_t i) const {
        return static_cast<StructData<1,Vector<T2,3>>&>(this->data[i]).data;
    }
};

#endif // VMD_H
