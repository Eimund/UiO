
/*
 *  FYS4150 - Computational Physics - Project 5
 *
 *  Written by:         Eimund Smestad
 *  Candidate number:   68
 *
 *  14.11.2014
 *
 *  c++11 compiler
 */

#include <cmath>
#include "time.h"
#include <chrono>
#include <fstream>
#include <random>

#include "Array.h"

using namespace std;

#define FLOAT double

/*template<class T> struct Chain {
    unsigned int N;
    T element;
    Chain<T>* owner;
    Chain<T>* prev;
    Chain<T>* next;
    Chain() {
        N = 0;
        owner = this;
        prev = this;
        next = this;
    }
    ~Chain() {
        if(N) {
            Chain<T>* p = this;
            for(int i = 1; i < N; i++)
                p = p->next;
            for(int i = 0; i < N; i++) {
                delete p->next;
                p = p->prev;
            }
        }
    }
    void Add(T element) {       // Add next
        if(next == this) {
            next = new Chain<T>;
            next->next = next;
        }
        else {
            next->prev = new Chain<T>;
            next->prev->next = next;
            next = next->prev;
        }
        owner->N++;
        next->owner = owner;
        next->prev = this;
        next->element = element;
    }
    void Remove() {             // Remove next
        if(next->next == next) {
            delete next;
            next = this;
        } else {
            next = next->next;
            delete next->prev;
            next->prev = this;
        }
        owner->N--;
    }
};

template<typename T, unsigned int D> struct Pointer : Pointer<T*,D-1> {
};
template<typename T> struct Pointer<T,0> {
    typedef T Type;
};

template<typename T, unsigned int D> struct Space {
    static typename Pointer<T,D>::Type Allocate(unsigned int n[D]) {
        typename Pointer<T,D>::Type array = new typename Pointer<T,D-1>::Type[n[0]];
        for(int i = 0; i < n[0]; i++)
            array[i] = Space<T,D-1>::Allocate(&n[1]);
        return array;
    }
    static void Deallocate(typename Pointer<T,D>::Type array, unsigned int n[D]) {
        for(int i = 0; i < n[0]; i++)
            Space<T,D-1>::Deallocate(array[i], &n[1]);
        delete [] array;
    }
    static void DeRange(T** array) {
        for(int i = 0; i < D; i++)
            delete [] array[i];
        delete [] array;
    }
    static T** Range(T lower[D], T upper[D], unsigned int n[D]) {
        T** array = new T*[D];
        for(int i = 0; i < D; i++) {
            T step = (upper[i]-lower[i])/(n[i]-1);
            array[i] = new T[n[i]];
            array[i][0] = lower[i];
            for(int j = 1; j < n[i]; j++)
                array[i][j] = array[i][j-1] + step;
        }
        return array;
    }
    static void ToFile(ofstream& file, typename Pointer<T,D>::Type array, unsigned int n[D]) {
        for(int i = 0; i < n[0]; i++) {
            if(D == 1 && i)
                file << '\t';
            else if(D > 1 && i)
                file << endl;

            Space<T,D-1>::ToFile(file, array[i], &n[1]);
        }
    }
};
template<typename T> struct Space<T,0> {
    static T Allocate(unsigned int[0]) {
        return T();
    }
    static void Deallocate(T, unsigned int[0]) {
    }
    static void ToFile(ofstream& file, T array, unsigned int[0]) {
        file << array;
    }
};
template<typename T> struct Space<T*,0> {
    static T* Allocate(unsigned int[0]) {
        return new T;
    }
    static void Deallocate(T* array, unsigned int[0]) {
        delete array;
    }
    static void ToFile(ofstream& file, T* array, unsigned int[0]) {
        file << *array;
    }
};

template<typename T, unsigned int D> struct Vector {
    T element[D];
    operator T*() {
        return element;
    }
    Vector & operator= (const T other[D]) {
        for(int i = 0; i < D; i++)
            element[i] = other[i];
    }
}; */

template<typename T> Chain<Vector<T,1>> Diffusion1D_MonteCarlo(unsigned int N, T dt, T t, T d, T (*step)());
template<typename T> Vector<T,2> Diffusion2D_alpha(T alpha, T d[2], unsigned int n[2]);
template<typename T> T Diffusion2D_deltaT(T alpha, T dx2[2]);
template<typename T> Vector<T,2> Diffusion2D_deltaX2(T d[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Exact(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]);
template<typename T> T** Diffusion2D_Explicit(T alpha, T t, T d[2], T w[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Initialize(unsigned int m[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Jacobi(T theta, T alpha, T t,T d[2], T w[2], unsigned int n[2]);
template<typename T> Vector<unsigned int,2> Diffusion2D_Source(T d, T w[2], unsigned int n);
template<typename T> T StepUniform();

int main() {

    ofstream file;

    FLOAT u0 = 1;
    FLOAT D = 1;
    FLOAT t = 0.1;
    Vector<unsigned int,2> n = {100,100};
    Vector<unsigned int,2> ne = {100,100};
    Vector<FLOAT,2> d = {1,10};
    Vector<FLOAT,2> w = {7,8};
    auto x = Space<FLOAT,2>::Range(ARRAYLIST(FLOAT,0,0),d,n);

    file.open("Exact_2D.dat");
    Space<FLOAT,1>::ToFile(file, x[0], n);
    file << endl;
    Space<FLOAT,1>::ToFile(file, x[1], &n[1]);
    file << endl;
    auto u = Diffusion2D_Exact<FLOAT>(u0, D, t, d, w, x, n, ne);
    Space<FLOAT,2>::ToFile(file, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("Explicit_2D.dat");
    Space<FLOAT,1>::ToFile(file, x[0], n);
    file << endl;
    Space<FLOAT,1>::ToFile(file, x[1], &n[1]);
    file << endl;
    u = Diffusion2D_Explicit<FLOAT>(0.1, t, d, w, n);
    Space<FLOAT,2>::ToFile(file, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("Jacobi_2D.dat");
    Space<FLOAT,1>::ToFile(file, x[0], n);
    file << endl;
    Space<FLOAT,1>::ToFile(file, x[1], &n[1]);
    file << endl;
    u = Diffusion2D_Jacobi<FLOAT>(1, 0.1, t, d, w, n);
    Space<FLOAT,2>::ToFile(file, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    auto particles = Diffusion1D_MonteCarlo<FLOAT>(10, 1e-2, t, d[0], &StepUniform);

    Space<FLOAT,2>::DeRange(x);

    return 0;
}

template<typename T> Chain<Vector<T,1>> Diffusion1D_MonteCarlo(unsigned int N, T dt, T t, T d, T (*step)()) {

    Chain<Vector<T,1>> particles, *p;
    Vector<T,1> vec = {0};
    for(int i = 0; i < N; i++)
        particles.Add(vec);                                 // Initial particles

    T val;
    unsigned int n = t / dt;
    T dx = sqrt(2*dt);
    uniform_real_distribution<T> pdf;                       // Distribution for accepting a move
    default_random_engine rng;                              // Set RNG seeding value
    rng.seed(chrono::high_resolution_clock::now().time_since_epoch().count());

    for(int i = 0; i < n; i++) {                            // Loop of timesteps
        p = &particles;
        for(int j = 0; j < particles.N; j++) {              // Loop of particles
            val = pdf(rng);                                 // Random number
            p = p->next;                                    // Next particle

            if(val <= 0.5) {                                // Sampling rule

                val = p->element.element[0] - dx * step();  // Calculate move
                if(p->element.element[0] >= 0)              // Valid move
                    p->element.element[0] = val;

            } else {

                if(p->element.element[0] == 0) {            // Particle at the source
                    p->prev->Add(p->element);               // Add new particle at prev
                    j++;
                }
                p->element.element[0] += dx * step();       // Move particle
                if(p->element.element[0] >= d) {            // Remove particle
                    p = p->prev;                            // into the postsynaptic
                    p->Remove();
                    j--;
                }
            }
        }
    }

    return particles;
}

template<typename T> Vector<T,2> Diffusion2D_alpha(T dt, T dx2[2]) {
    Vector<T,2> vec;
    vec.element[0] = dt/dx2[0];
    vec.element[1] = dt/dx2[1];
    return vec;
}

template<typename T> T Diffusion2D_deltaT(T alpha, T dx2[2]) {
    T dt0 = alpha * dx2[0];
    T dt1 = alpha * dx2[1];

    if(dt0 < dt1)
        return dt0;
    return dt1;
}

template<typename T> Vector<T,2> Diffusion2D_deltaX2(T d[2], unsigned int n[2]) {
    Vector<T,2> x2;
    x2.element[0] = d[0]/n[0];
    x2.element[0] *= x2.element[0];
    x2.element[1] = d[1]/n[1];
    x2.element[1] *= x2.element[1];
    return x2;
}

template<typename T> T** Diffusion2D_Exact(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]) {
    T omega, a, f;
    T t1 = - D*t/(w[0]*w[0]);
    T t2 = - D*t/((w[1]-d[1])*(w[1]-d[1]));
    T t3 = - D*t/(d[0]*d[0]);
    auto u = Space<T,2>::Allocate(n);

    for(int i = 0; i < n[0]; i++) {
        for(int j = 0; j < n[1]; j++) {

            if(x[1][j] < w[0]) {
                a = 0;
                f = 1-x[1][j]/w[0];
                for(int k=1; k < ne[1]; k++) {
                    omega = M_PI* k;
                    a -= (T)2/omega * sin(omega*f) * exp(omega*omega*t1);
                }
                a += 1 - f;
            }
            else if(x[1][j] > w[1]) {
                a = (x[1][j]-d[1])/(w[1]-d[1]);
                f = (w[1]-x[1][j])/(w[1]-d[1]);
                for(int k=1; k < ne[1]; k++) {
                    omega = M_PI* k;
                    a -= (T)2/omega * sin(omega*f) * exp(omega*omega*t2);
                }
            }
            else
                a = 1;

            u[i][j] = 0;
            f = x[0][i]/d[0];
            for(int k = 1; k < ne[0]; k++) {
                omega = M_PI* k;
                u[i][j] -= (T)2/omega * sin(omega*f) * exp(omega*omega*t3);
            }
            u[i][j] = u0 * a * (1 - f + u[i][j]);
        }
    }
    return u;
}

template<typename T> T** Diffusion2D_Explicit(T alpha, T t, T d[2], T w[2], unsigned int n[2]) {
    Vector<unsigned int,2> m = Diffusion2D_Source(d[1], w, n[1]);
    T** u = Diffusion2D_Initialize<T>(m, n);
    T** v = Diffusion2D_Initialize<T>(m, n);
    T** z;
    Vector<T,2> dx2 = Diffusion2D_deltaX2(d, n);
    T dt = Diffusion2D_deltaT<T>(alpha, dx2);
    Vector<T,2> a = Diffusion2D_alpha<T>(dt, dx2);
    unsigned int nt = t/dt;

    n[0]--;
    n[1]--;
    for(int i = 0, j, k; i < nt; i++) {

        for(k = 1; k < m.element[0]; k++) {
            v[0][k] = u[0][k];
            v[0][k] += a.element[1]*(u[0][k+1] - 2.0*u[0][k] + u[0][k-1]);
        }
        for(k = m.element[1]+1; k < n[1]; k++) {
            v[0][k] = u[0][k];
            v[0][k] += a.element[1]*(u[0][k+1] - 2.0*u[0][k] + u[0][k-1]);
        }

        for(j = 1; j < n[0]; j++) {
            for(k = 1; k < n[1]; k++) {
                v[j][k] = u[j][k];
                v[j][k] += a.element[0]*(u[j+1][k] - 2.0*u[j][k] + u[j-1][k]);
                v[j][k] += a.element[1]*(u[j][k+1] - 2.0*u[j][k] + u[j][k-1]);
            }
        }
        z = u;
        u = v;
        v = z;
    }
    n[0]++;
    n[1]++;
    Space<T,2>::Deallocate(v, n);

    return u;
}

template<typename T> T** Diffusion2D_Initialize(unsigned int m[2], unsigned int n[2]) {
    auto u = Space<T,2>::Allocate(n);
    for(int i = 0; i < n[0]; i++) {
        for(int j = 0; j < n[1]; j++)
            u[i][j] = 0;
    }
    for(int i = m[0]; i <= m[1]; i++)
        u[0][i] = 1;

    return u;
}

template<typename T> T** Diffusion2D_Jacobi(T theta, T alpha, T t,T d[2], T w[2], unsigned int n[2]) {
    Vector<unsigned int,2> m = Diffusion2D_Source(d[1], w, n[1]);
    T** u = Diffusion2D_Initialize<T>(m, n);
    T** v = Diffusion2D_Initialize<T>(m, n);
    auto c = Space<T,2>::Allocate(n);
    T** z;
    Vector<T,2> dx2 = Diffusion2D_deltaX2(d, n);
    T dt = Diffusion2D_deltaT<T>(alpha, dx2);
    Vector<T,2> a = Diffusion2D_alpha<T>(dt, dx2);
    unsigned int nt = t/dt;

    n[0]--;
    n[1]--;

    T diff;
    T c0 = theta * a.element[1];
    T c1 = 1 + 2 * c0;
    T c2 = theta * a.element[0];
    T c3 = 1 + 2 * (c0 + c2);
    T c4 = a.element[1] - c0;
    T c5 = 1 - 2 * c4;
    T c6 = a.element[0] - c2;
    T c7 = (c5 - 2 * c6) / c3;
    T c8 = c6 / c3;
    T c9 = c4 / c3;
    c5 = c5 / c1;
    c6 = c4 / c1;
    c1 = c0 / c1;
    c4 = c0 / c3;
    c3 = c2 / c3;


    for(int i = 0, j, k; i < nt; i++) {

        for(k = 1; k < m.element[0]; k++)
            c[0][k] = c5 * u[0][k] + c6 * (u[0][k+1] + u[0][k-1]);
        for(k = m.element[1]+1; k < n[1]; k++)
            c[0][k] = c5 * u[0][k] + c6 * (u[0][k+1] + u[0][k-1]);
        for(j = 1; j < n[0]; j++) {
            for(k = 1; k < n[1]; k++)
                c[j][k] = c7 * u[j][k] + c8 * (u[j+1][k] + u[j-1][k]) + c9 * (u[j][k+1] + u[j][k-1]);
        }

        do {
            diff = 0;

            for(k = 1; k < m.element[0]; k++) {
                v[0][k] = c[0][k] + c1 * (u[0][k+1] + u[0][k-1]);
                diff += abs(v[0][k] - u[0][k]);
            }
            for(k = m.element[1]+1; k < n[1]; k++) {
                v[0][k] = c[0][k] + c1 * (u[0][k+1] + u[0][k-1]);
                diff += abs(v[0][k] - u[0][k]);
            }
            for(j = 1; j < n[0]; j++) {
                for(k = 1; k < n[1]; k++) {
                    v[j][k] = c[j][k] + c3 * (u[j+1][k] + u[j-1][k]) + c4 * (u[j+1][k] + u[j-1][k]);
                    diff += abs(v[j][k] - u[j][k]);
                }
            }

            diff /= n[0] * n[1];

            z = u;
            u = v;
            v = z;

        } while(diff > 0.00001);
    }

    n[0]++;
    n[1]++;
    Space<T,2>::Deallocate(v, n);
    Space<T,2>::Deallocate(c, n);

    return u;
}

template<typename T> Vector<unsigned int,2> Diffusion2D_Source(T d, T w[2], unsigned int n) {
    Vector<unsigned int,2> m;
    m.element[0] = w[0]/d *(n-1) + (T)0.5;
    m.element[1] = w[1]/d *(n-1) + (T)0.5;
    return m;
}

template<typename T> T StepUniform() {
        return (T)1;
}
