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
        public: inline static U DT(U dt) {
            return dt/2;
        }
        public: inline static U Step(unsigned int i, U* (*f)(U, U*, U*), U** x0, U** v0, U t, U dt, unsigned int dim) {
            U* x = new T[dim];
            U* v = new T[dim];
            U* x2 = new T[dim];
            U* v2 = new T[dim];
            unsigned int j;
            for(j = 0; j < dim; j++) {
                x[j] = x0[j][i-1];
                v[j] = v0[j][i-1];
                x2[j] = x[j];
                v2[j] = v[j];
            }
            x2 = f(t, x2, v2);
            for(j = 0; j < dim; j++) {
                x0[j][i] = x2[j];
                v0[j][i] = v2[j];
                x2[j] = x[j] + x2[j]*dt;
                v2[j] = v[j] + v2[j]*dt;
            }
            t += dt;
            x2 = f(t, x2, v2);
            for(j = 0; j < dim; j++) {
                x0[j][i] += 2*x2[j];
                v0[j][i] += 2*v2[j];
                x2[j] = x[j] + x2[j]*dt;
                v2[j] = v[j] + v2[j]*dt;
            }
            x2 = f(t, x2, v2);
            t += dt;
            dt += dt;
            for(j = 0; j < dim; j++) {
                x0[j][i] += 2*x2[j];
                v0[j][i] += 2*v2[j];
                x2[j] = x[j] + x2[j]*dt;
                v2[j] = v[j] + v2[j]*dt;
            }
            x2 = f(t, x2, v2);
            for(j = 0; j < dim; j++) {
                x0[j][i] += x2[j];
                x0[j][i] *= dt/6;
                x0[j][i] += x0[j][i-1];
                v0[j][i] += x2[j];
                v0[j][i] *= dt/6;
                v0[j][i] += v0[j][i-1];
            }
            delete [] x;
            delete [] v;
            delete [] x2;
            delete [] v2;
            return t;
        }
        public: static U* func(U* (*f)(U, U*), U t, U* x) {
            return x;
        }
    };
    private: template<class U> class Diff<DifferentialType::Verlet, U> {
        public: inline static U DT(U dt) {
            return dt;
        }
        public: inline static U Step(unsigned int i, U* (*f)(U, U*), U** x0, U t, U dt, unsigned int dim) {
            U* x2 = new U[dim];
            for(unsigned int j = 0; j < dim; j++)
                x2[j] = x0[j][i-1];
            x2 = f(t, x2);
            for(unsigned int j = 0; j < dim; j++)
                x0[j][i] = 2*x0[j][i-1]-x0[j][i-2]+dt*dt*x2[j];
            delete [] x2;
            return t+dt;
        }
        public: static U* func(U* (*f)(U, U*), U t, U* x) {
            return f(t,x);
        }
    };
    public: template<DifferentialType Type> T* func(T t, T* x) {
        return Diff<Type,T>::func(f, t, x);
    }
    public: template<DifferentialType Type> T Step(unsigned int i, T* (*f)(T, T*), T** x0, T t, T dt, unsigned int dim) {
        this->f = f;
        return Diff<Type,T>::Step(i, &Differential_2<T>::func<Type>, x0, t, dt, dim);
    }
    public: template<DifferentialType Type> T* Solve(T* (*f)(T, T*), T** x0, T dt, unsigned int dim, unsigned int step) {
        T* t = new T[step];
        t[0] = -dt;
        t[1] = 0;
        dt = Diff<Type, T>::DT(dt);
        this->f = f;
        for(unsigned int i = 2; i < step; i++)
            t[i] = Diff<Type,T>::Step(i, &Differential_2<T>::func<Type>, x0, t[i-1], dt, dim);
        return t;
    }
};

#endif // DIFFERENTIAL_H
