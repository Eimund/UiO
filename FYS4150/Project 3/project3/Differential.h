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

template<class T>  class Differential_2 {
    private: T* (*f)(T, T*);
    private: template<DifferentialType Type, class C> class Diff;
    private: template<class C> class Diff<DifferentialType::RK4, C> {
    public: inline static void Solve(C* owner, T* (C::*f)(T*), T* t, T** x0, T** v0, T dt, unsigned int dim, unsigned int step) {
            T* a;
            T* x = new T[dim];
            T* x2 = new T[dim];
            T* v = new T[dim];
            T* v2 = new T[dim];
            T dt2 = dt/2;
            T dt6 = dt/6;
            for(unsigned int i = 1, i1 = 0, j; i < step; i++, i1++) {
                for(j = 0; j < dim; j++) {
                    x[j] = x0[j][i1];
                    v[j] = v0[j][i1];
                }
                a = (owner->*f)(x);
                for(j = 0; j < dim; j++) {
                    x0[j][i] = v[j];
                    v0[j][i] = a[j];
                    x2[j] = x[j] + v[j]*dt2;
                    v2[j] = v[j] + a[j]*dt2;
                }
                a = (owner->*f)(x2);
                for(j = 0; j < dim; j++) {
                    x0[j][i] += 2*v2[j];
                    v0[j][i] += 2*a[j];
                    x2[j] = x[j] + v2[j]*dt2;
                    v2[j] = v[j] + a[j]*dt2;
                }
                a = (owner->*f)(x2);
                for(j = 0; j < dim; j++) {
                    x0[j][i] += 2*v2[j];
                    v0[j][i] += 2*a[j];
                    x2[j] = x[j] + v2[j]*dt;
                    v2[j] = v[j] + a[j]*dt;
                }
                a = (owner->*f)(x2);
                for(j = 0; j < dim; j++) {
                    x0[j][i] += v2[j];
                    x0[j][i] *= dt6;
                    x0[j][i] += x[j];
                    v0[j][i] += a[j];
                    v0[j][i] *= dt6;
                    v0[j][i] += v[j];
                }
                t[i] = t[i1]+dt;
            }
            delete [] x;
            delete [] v;
            delete [] x2;
            delete [] v2;
        }
    };
    private: template<class C> class Diff<DifferentialType::Verlet, C> {
        public: inline static void Solve(C* owner, T* (C::*f)(T*), T* t, T** x0, T** v0, T dt, unsigned int dim, unsigned int step) {
            T* a;
            T* x2 = new T[dim];

            // Initial step
            for(unsigned int j = 0; j < dim; j++)
                x2[j] = x0[j][0];
            a = (owner->*f)(x2);
            for(unsigned int j = 0; j < dim; j++)
                x0[j][1] = x2[j] + dt*(v0[j][0] + dt*a[j]);

            // Continuing step
            for(unsigned int i = 2, i1=1, i2=0, j; i < step; i++, i1++, i2++) {
                for(j = 0; j < dim; j++)
                    x2[j] = x0[j][i1];
                a = (owner->*f)(x2);
                for(j = 0; j < dim; j++) {
                    x0[j][i] = 2*x2[j] - x0[j][i2] + dt*dt*a[j];
                    v0[j][i1] = (x0[j][i]-x0[j][i2])/(2*dt);
                }
                t[i] = t[i1]+dt;
            }
            unsigned int i1 = step-1;
            unsigned int i2 = step-2;
            for(unsigned int j = 0; j < dim; j++)
                v0[j][i1] = v0[j][i2] + a[j]*dt;
            delete [] x2;
        }
    };
    public: template<DifferentialType Type, class C> void Solve(C* owner, T* (C::*f)(T*), T* t, T** x0, T** v0, T dt, unsigned int dim, unsigned int step) {
        t[0] = 0;
        t[1] = dt;
        Diff<Type,C>::Solve(owner, f, t, x0, v0, dt, dim, step);
    }
};

#endif // DIFFERENTIAL_H
