/*
* FYS4150 - Computational Physics - Project 3
*
* Written by: Eimund Smestad
*
* 15.10.2014
*
* c11 compiler
*/

#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H

enum class DifferentialType {
    RK4,
    Verlet
};

class Differential_1 {
    private: template<DifferentialType Type, class T> class Diff;
    private: template<class T> class Diff<DifferentialType::RK4, T> {
        public: inline static T Step(unsigned int i, T** x0, T* (*f)(T, T*), T t, T dt, unsigned int dim) {
            T* x = new T[dim];
            T* x2 = new T[dim];
            unsigned int j;
            for(j = 0; j < dim; j++) {
                x[j] = x0[j][i-1];
                x2[j] = x[j];
            }
            x2 = f(t, x2);
            for(j = 0; j < dim; j++) {
                x0[j][i] = x2[j];
                x2[j] = x[j] + x2[j]*dt;
            }
            t += dt;
            x2 = f(t, x2);
            for(j = 0; j < dim; j++) {
                x0[j][i] += 2*x2[j];
                x2[j] = x[j] + x2[j]*dt;
            }
            x2 = f(t, x2);
            t += dt;
            dt += dt;
            for(j = 0; j < dim; j++) {
                x0[j][i] += 2*x2[j];
                x2[j] = x[j] + x2[j]*dt;
            }
            x2 = f(t, x2);
            for(j = 0; j < dim; j++) {
                x0[j][i] += x2[j];
                x0[j][i] *= dt/6;
                x0[j][i] += x0[j][i-1];
            }
            delete [] x;
            delete [] x2;
            return t;
        }
    };
    public: template<DifferentialType Type, class T> inline T Step(unsigned int i, T** x0, T* (*f)(T, T*), T t, T dt, unsigned int dim) {
        return Diff<Type,T>::Step(i, x0, f, t, dt, dim);
    }
};
class Differential_2 : public Differential_1 {
    private: template<DifferentialType Type, class T> class Diff;
    private: template<class T> class Diff<DifferentialType::RK4, T> {
        public: inline static T DT(T dt) {
            return dt/2;
        }
        public: inline static void Step(unsigned int i, T** x0, T* (*f)(T, T*), T t, T dt, unsigned int dim) {
        }
    };
    private: template<class T> class Diff<DifferentialType::Verlet, T> {
        public: inline static T DT(T dt) {
            return dt;
        }
        public: inline static T Step(unsigned int i, T** x0, T* (*f)(T, T*), T t, T dt, unsigned int dim) {
            T* x2 = new T[dim];
            for(unsigned int j = 0; j < dim; j++)
                x2[j] = x0[j][i-1];
            x2 = f(t, x2);
            for(unsigned int j = 0; j < dim; j++)
                x0[j][i] = 2*x0[j][i-1]-x0[j][i-2]+dt*dt*x2[j];
            delete [] x2;
            return t+dt;
        }
    };
    public: template<DifferentialType Type, class T> T Step(unsigned int i, T** x0, T* (*f)(T, T*), T t, T dt, unsigned int dim) {
        return Diff<Type,T>::Step(i, x0, f, t, dt, dim);
    }
    public: template<DifferentialType Type, class T> T* Solve(T** x0, T* (*f)(T, T*), T dt, unsigned int dim, unsigned int step) {
        T* t = new T[step];
        t[0] = 0;
        t[1] = dt;
        dt = Diff<Type, T>::DT(dt);
        for(unsigned int i = 2; i < step; i++)
            t[i] = Diff<Type, T>::Step(i, x0, f, t[i-1], dt, dim);
        return t;
    }
};

#endif // DIFFERENTIAL_H
