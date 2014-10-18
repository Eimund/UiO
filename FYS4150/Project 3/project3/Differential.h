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

template<class T> class Differential_1 {
    private: T* (*f)(T, T*);
    private: template<DifferentialType Type, class U> struct Diff;
    private: template<class U> struct Diff<DifferentialType::RK4, U> {
        public: inline static U Step( unsigned int i, U* (*f)(U, U*), U** x0, U t, U dt, unsigned int dim) {
            U* x = new T[dim];
            U* x2 = new T[dim];
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
        public: static U* func(U* (*f)(U, U*), U t, U* x) {
            return f(t,x);
        }
    };
    public: template<DifferentialType Type> T* func(T t, T* x) {
        return Diff<Type,T>::func(f, t, x);
    }
    public: template<DifferentialType Type> inline T Step(unsigned int i, T* (*f)(T, T*), T** x0, T t, T dt, unsigned int dim) {
        this->f = f;
        return Diff<Type,T>::Step(i, &Differential_1<T>::func<Type>, x0, t, dt, dim);
    }
};
template<class T>  class Differential_2 {
    private: T* (*f)(T, T*);
    private: template<DifferentialType Type, class U> class Diff;
    private: template<class U> class Diff<DifferentialType::RK4, U> {
        public: inline static void Solve(U* t, T* (*f)(T*), T** x0, T dt, unsigned int dim, unsigned int step) {
            U* a;
            U* x = new T[dim];
            U* x2 = new T[dim];
            U* v = new T[dim];
            for(unsigned int i = 0; i < dim; i++)
                v[i] = (x[i][1] - x[i][0])/dt;
            U* v1 = new T[dim];
            U* v2 = new T[dim];
            U dt2 = dt/2;
            U dt6 = dt/6;
            for(unsigned int i = 2, i1 = 1, j; i < step; i++, i1++) {
                for(j = 0; j < dim; j++)
                    x[j] = x[j][i1];
                a = f(x);
                for(j = 0; j < dim; j++) {
                    x[j][i] = v[j];
                    v1[i] = a[j];
                    x2[j] = x[j] + v[j]*dt2;
                    v2[j] = v[j] + a[j]*dt2;
                }
                a = f(x2);
                for(j = 0; j < dim; j++) {
                    x[j][i] += 2*v2[j];
                    v1[i] += 2*a[j];
                    x2[j] = x[j] + v2[j]*dt2;
                    v2[j] = v[j] + a[j]*dt2;
                }
                a = f(x2);
                for(j = 0; j < dim; j++) {
                    x[j][i] += 2*v2[j];
                    v1[i] += 2*a[j];
                    x2[j] = x[j] + v2[j]*dt;
                    v2[j] = v[j] + a[j]*dt;
                }
                a = f(x2);
                for(j = 0; j < dim; j++) {
                    x[j][i] += v2[j];
                    x[j][i] *= dt6;
                    x[j][i] += x[j];
                    v1[j] += a[j];
                    v1[j] *= dt6;
                    v[j] += v1[j];
                }
                t[i] = t[i1]+dt;
            }
            delete [] x;
            delete [] x2;
            delete [] v;
            delete [] v1;
            delete [] v2;
        }
    };
    private: template<class U> class Diff<DifferentialType::Verlet, U> {
        public: inline static void Solve(U* t, T* (*f)(T*), T** x0, T dt, unsigned int dim, unsigned int step) {
            U* a;
            U* x2 = new U[dim];
            for(unsigned int i = 2, i1=1, i2=0, j; i < step; i++, i1++, i2++) {
                for(j = 0; j < dim; j++)
                    x2[j] = x0[j][i1];
                a = f(x2);
                for(j = 0; j < dim; j++)
                    x0[j][i] = 2*x2[j] - x0[j][i2] + dt*dt*a[j];
                t[i] = t[i1]+dt;
            }
            delete [] x2;
        }
    };
    public: template<DifferentialType Type> T* Solve(T* (*f)(T*), T** x0, T dt, unsigned int dim, unsigned int step) {
        T* t = new T[step];
        t[0] = -dt;
        t[1] = 0;
        Diff<Type,T>::Solve(t, f, x0, dt, dim, step);
        return t;
    }
};

#endif // DIFFERENTIAL_H
