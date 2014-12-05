
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
#include <cstdarg>

#include "Array.h"
#include "Delegate.h"
#include "Experiment.h"

using namespace std;

#define FLOAT double

template<typename T> struct Step {
    uniform_real_distribution<T> pdf;
    default_random_engine rng;
    Step() : pdf(0.0,1.0) {
        rng.seed(chrono::high_resolution_clock::now().time_since_epoch().count());
    }
    T Uniform() {
        return (T)1;
    }
    T Gaussian() {
        return sqrt(-log(pdf(rng)));
    }
};

template<typename C, typename T> Chain<Vector<T,1>> Diffusion1D_MonteCarlo(T t, T dt, T d, Delegate<C,T> step);
template<typename T> T* Diffusion1D_Metropolis(unsigned long N, unsigned int n);
template<typename T> Vector<T,2> Diffusion2D_alpha(T alpha, T d[2], unsigned int n[2]);
template<typename T> T Diffusion2D_deltaT(T alpha, T dx2[2]);
template<typename T> Vector<T,2> Diffusion2D_deltaX2(T d[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Exact(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]);
template<typename T> T** Diffusion2D_Exact2(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]);
template<typename T> T** Diffusion2D_Explicit(T alpha, T t, T d[2], T w[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Explicit2(T alpha, T t, T d[2], T w[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Initialize(unsigned int m[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Jacobi(T theta, T alpha, T t,T d[2], T w[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Jacobi2(T theta, T alpha, T t,T d[2], T w[2], unsigned int n[2]);
template<typename C, typename T> Chain<Vector<T,2>> Diffusion2D_MonteCarlo(T t, T dt, T d[2], T w[2], Delegate<C,T> step);
template<typename T> Vector<unsigned int,2> Diffusion2D_Source(T d, T w[2], unsigned int n);
template<typename T> T StepUniform();
template<typename T> T StepGaussian();

int main() {

    ofstream file;

    FLOAT u0 = 1;
    FLOAT D = 1;
    FLOAT t = 1;
    Vector<unsigned int,2> n = {100,100};
    Vector<unsigned int,2> ne = {100,100};
    Vector<FLOAT,2> d = {1,10};
    Vector<FLOAT,2> w = {7,8};
    auto x = Space<FLOAT,2>::Range(ARRAYLIST(FLOAT,0,0), d, n);

    file.open("Exact_2D.dat");
    auto u = Diffusion2D_Exact<FLOAT>(u0, D, t, d, w, x, n, ne);
    Space<FLOAT,2>::ToFile(file, x, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("Explicit_2D.dat");
    u = Diffusion2D_Explicit2<FLOAT>(0.1, t, d, w, n);
    Space<FLOAT,2>::ToFile(file, x, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("Jacobi_2D.dat");
    u = Diffusion2D_Jacobi<FLOAT>(1, 0.1, t, d, w, n);
    Space<FLOAT,2>::ToFile(file, x, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("MonteCarlo_1D_Uniform.dat");
    Step<FLOAT> step;
    Delegate<Step<FLOAT>,FLOAT> step_uniform(&step, &Step<FLOAT>::Uniform);
    Delegate<void, Chain<Vector<FLOAT,1>>, FLOAT, FLOAT, FLOAT, decltype(step_uniform)> Diff1D_MC(&Diffusion1D_MonteCarlo<Step<FLOAT>,FLOAT>);
    auto space = Experiment(1000u, x, n, Diff1D_MC, t, 1e-3, d[0], step_uniform);
    Space<FLOAT,1>::ToFile(file, x, space, n);
    Space<FLOAT,1>::Deallocate(space, n);
    file.close();

    file.open("MonteCarlo_1D_Gaussian.dat");
    Delegate<Step<FLOAT>,FLOAT> step_gaussian(&step, &Step<FLOAT>::Gaussian);
    space = Experiment(1000u, x, n, Diff1D_MC, t, 1e-3, d[0], step_gaussian);
    Space<FLOAT,1>::ToFile(file, x, space, n);
    Space<FLOAT,1>::Deallocate(space, n);
    file.close();

    file.open("MonteCarlo_2D_Gaussian.dat");
    Delegate<void, Chain<Vector<FLOAT,2>>, FLOAT, FLOAT, FLOAT*, FLOAT*, decltype(step_uniform)> Diff2D_MC(&Diffusion2D_MonteCarlo<Step<FLOAT>,FLOAT>);
    auto space2d = Experiment(1000u, x, n, Diff2D_MC, t, 1e-3, d.e, w.e, step_gaussian);
    Space<FLOAT,2>::ToFile(file, x, space2d, n);
    Space<FLOAT,2>::Deallocate(space2d, n);
    file.close();

    file.open("Explicit_2D_2.dat");
    u = Diffusion2D_Explicit<FLOAT>(0.1, t, d, w, n);
    Space<FLOAT,2>::ToFile(file, x, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("Jacobi_2D_2.dat");
    u = Diffusion2D_Jacobi2<FLOAT>(1, 0.1, t, d, w, n);
    Space<FLOAT,2>::ToFile(file, x, u, n);
    Space<FLOAT,2>::Deallocate(u, n);
    file.close();

    file.open("Metropolis_1D.dat");
    space = Diffusion1D_Metropolis<FLOAT>(10000, n[0]);
    Space<FLOAT,1>::ToFile(file, x, space, n);
    Space<FLOAT,1>::Deallocate(space, n);
    file.close();

    Space<FLOAT,2>::DeRange(x);

    int i = 0;

    return 0;
}

template<typename C, typename T> Chain<Vector<T,1>> Diffusion1D_MonteCarlo(T t, T dt, T d, Delegate<C,T> step) {
    T val;
    unsigned int n = t / dt;
    T dx = sqrt(2*dt);
    uniform_real_distribution<T> pdf(0.0,1.0);      // Distribution for accepting a move
    default_random_engine rng;                      // Set RNG seeding value
    rng.seed(chrono::high_resolution_clock::now().time_since_epoch().count());

    Chain<Vector<T,1>> particles, *p;
    Vector<T,1> vec = {0};
    particles.Add(vec);                             // Seeding particle

    for(int i = 0; i < n; i++) {                    // Loop of timesteps
        p = &particles;
        for(int j = 0; j < particles.N; j++) {      // Loop of particles
            val = pdf(rng);                         // Random number
            p = p->next;                            // Next particle

            if(val <= 0.5) {                        // Sampling rule

                val = p->e.e[0] - dx * step();      // Calculate backward move
                if(val > 0)                         // Valid move
                    p->e.e[0] = val;

            } else {

                if(p->e.e[0] == 0) {                // Particle at the source
                    p->prev->Add(p->e);             // Add new particle at prev
                    j++;
                }
                p->e.e[0] += dx * step();           // Forward move particle
                if(p->e.e[0] >= d) {                // Remove particle
                    p = p->prev;                    // into the postsynaptic
                    p->Remove();
                    j--;
                }
            }
        }
    }

    return particles;
}

template<typename T> T* Diffusion1D_Metropolis(unsigned long N, unsigned int n) {

    T val;
    auto s = Space<T,1>::Allocate(&n);
    s[0] = 1;

    unsigned int j = 0;
    for(unsigned long i = 0; i < N; i++) {
        for(j = 1; j < n; j++)
            Metropolis(0.5, s[j-1], 0.5, s[j]);
        s[0] = 1;
        s[n-1] = 0;
        for(j = n-1; j; j--)
            Metropolis(0.5, s[j], 0.5, s[j-1]);
        s[0] = 1;
        s[n-1] = 0;
    }
    return s;
}

template<typename T> Vector<T,2> Diffusion2D_alpha(T dt, T dx2[2]) {
    Vector<T,2> vec;
    vec.e[0] = dt/dx2[0];
    vec.e[1] = dt/dx2[1];
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
    x2.e[0] = d[0]/n[0];
    x2.e[0] *= x2.e[0];
    x2.e[1] = d[1]/n[1];
    x2.e[1] *= x2.e[1];
    return x2;
}

template<typename T> T** Diffusion2D_Exact(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]) {
    T omega, a, b, f;
    T t1 = - D*t/(3*w[0]*w[0]);
    T t2 = - D*t/(3*(w[1]-d[1])*(w[1]-d[1]));
    T t3 = 5*t1;
    T t4 = 5*t2;
    T t5 = - D*t/(d[0]*d[0]);
    auto u = Space<T,2>::Allocate(n);

    for(int i = 0; i < n[0]; i++) {
        for(int j = 0; j < n[1]; j++) {

            if(x[1][j] < w[0]) {
                a = 0;
                b = 0;
                f = 1-x[1][j]/w[0];
                for(int k=1; k < ne[1]; k++) {
                    omega = M_PI* k;
                    a -= (T)2/omega * sin(omega*f) * exp(omega*omega*t3);
                    b += (T)2/(sqrt(3)*omega) * exp(omega*omega*t1);
                }
                //b = 1;
                a += (1 - f) * b;
            }
            else if(x[1][j] > w[1]) {
                a = 0;
                b = 0;
                f = (w[1]-x[1][j])/(w[1]-d[1]);
                for(int k=1; k < ne[1]; k++) {
                    omega = M_PI* k;
                    a -= (T)2/omega * sin(omega*f) * exp(omega*omega*t4);
                    b += exp(omega*omega*t2);
                }
                b = 1;
                a += (x[1][j]-d[1])/(w[1]-d[1]) * b;
            }
            else
                a = 1;

            u[i][j] = 0;
            f = x[0][i]/d[0];
            for(int k = 1; k < ne[0]; k++) {
                omega = M_PI* k;
                u[i][j] -= (T)2/omega * sin(omega*f) * exp(omega*omega*t5);
            }
            u[i][j] = u0 * a * (1 - f + u[i][j]);
        }
    }
    return u;
}

template<typename T> T** Diffusion2D_Exact2(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]) {
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

        for(j = 1; j < n[0]; j++) {
            for(k = 1; k < n[1]; k++) {
                v[j][k] = u[j][k];
                v[j][k] += a.e[0]*(u[j+1][k] - 2.0*u[j][k] + u[j-1][k]);
                v[j][k] += a.e[1]*(u[j][k+1] - 2.0*u[j][k] + u[j][k-1]);
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
template<typename T> T** Diffusion2D_Explicit2(T alpha, T t, T d[2], T w[2], unsigned int n[2]) {
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

        for(k = 1; k < m.e[0]; k++) {
            v[0][k] = u[0][k];
            v[0][k] += a.e[1]*(u[0][k+1] - 2.0*u[0][k] + u[0][k-1]);
        }
        for(k = m.e[1]+1; k < n[1]; k++) {
            v[0][k] = u[0][k];
            v[0][k] += a.e[1]*(u[0][k+1] - 2.0*u[0][k] + u[0][k-1]);
        }

        for(j = 1; j < n[0]; j++) {
            for(k = 1; k < n[1]; k++) {
                v[j][k] = u[j][k];
                v[j][k] += a.e[0]*(u[j+1][k] - 2.0*u[j][k] + u[j-1][k]);
                v[j][k] += a.e[1]*(u[j][k+1] - 2.0*u[j][k] + u[j][k-1]);
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
    T c0 = theta * a.e[1];
    T c1 = 1 + 2 * c0;
    T c2 = theta * a.e[0];
    T c3 = 1 + 2 * (c0 + c2);
    T c4 = a.e[1] - c0;
    T c5 = 1 - 2 * c4;
    T c6 = a.e[0] - c2;
    T c7 = (c5 - 2 * c6) / c3;
    T c8 = c6 / c3;
    T c9 = c4 / c3;
    c5 = c5 / c1;
    c6 = c4 / c1;
    c1 = c0 / c1;
    c4 = c0 / c3;
    c3 = c2 / c3;


    for(int i = 0, j, k; i < nt; i++) {

        for(k = 1; k < m.e[0]; k++)
            c[0][k] = c5 * u[0][k] + c6 * (u[0][k+1] + u[0][k-1]);
        for(k = m.e[1]+1; k < n[1]; k++)
            c[0][k] = c5 * u[0][k] + c6 * (u[0][k+1] + u[0][k-1]);
        for(j = 1; j < n[0]; j++) {
            for(k = 1; k < n[1]; k++)
                c[j][k] = c7 * u[j][k] + c8 * (u[j+1][k] + u[j-1][k]) + c9 * (u[j][k+1] + u[j][k-1]);
        }

        do {
            diff = 0;

            for(k = 1; k < m.e[0]; k++) {
                v[0][k] = c[0][k] + c1 * (u[0][k+1] + u[0][k-1]);
                diff += abs(v[0][k] - u[0][k]);
            }
            for(k = m.e[1]+1; k < n[1]; k++) {
                v[0][k] = c[0][k] + c1 * (u[0][k+1] + u[0][k-1]);
                diff += abs(v[0][k] - u[0][k]);
            }
            for(j = 1; j < n[0]; j++) {
                for(k = 1; k < n[1]; k++) {
                    v[j][k] = c[j][k] + c3 * (u[j+1][k] + u[j-1][k]) + c4 * (u[j][k+1] + u[j][k-1]);
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

template<typename T> T** Diffusion2D_Jacobi2(T theta, T alpha, T t,T d[2], T w[2], unsigned int n[2]) {
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
    T c0 = theta * a.e[1];
    T c1 = 1 + 2 * c0;
    T c2 = theta * a.e[0];
    T c3 = 1 + 2 * (c0 + c2);
    T c4 = a.e[1] - c0;
    T c5 = 1 - 2 * c4;
    T c6 = a.e[0] - c2;
    T c7 = (c5 - 2 * c6) / c3;
    T c8 = c6 / c3;
    T c9 = c4 / c3;
    c5 = c5 / c1;
    c6 = c4 / c1;
    c1 = c0 / c1;
    c4 = c0 / c3;
    c3 = c2 / c3;


    for(int i = 0, j, k; i < nt; i++) {

        /*for(k = 1; k < m.e[0]; k++)
            c[0][k] = c5 * u[0][k] + c6 * (u[0][k+1] + u[0][k-1]);
        for(k = m.e[1]+1; k < n[1]; k++)
            c[0][k] = c5 * u[0][k] + c6 * (u[0][k+1] + u[0][k-1]);*/
        for(j = 1; j < n[0]; j++) {
            for(k = 1; k < n[1]; k++)
                c[j][k] = c7 * u[j][k] + c8 * (u[j+1][k] + u[j-1][k]) + c9 * (u[j][k+1] + u[j][k-1]);
        }

        do {
            diff = 0;

            /*for(k = 1; k < m.e[0]; k++) {
                v[0][k] = c[0][k] + c1 * (u[0][k+1] + u[0][k-1]);
                diff += abs(v[0][k] - u[0][k]);
            }
            for(k = m.e[1]+1; k < n[1]; k++) {
                v[0][k] = c[0][k] + c1 * (u[0][k+1] + u[0][k-1]);
                diff += abs(v[0][k] - u[0][k]);
            }*/
            for(j = 1; j < n[0]; j++) {
                for(k = 1; k < n[1]; k++) {
                    v[j][k] = c[j][k] + c3 * (u[j+1][k] + u[j-1][k]) + c4 * (u[j][k+1] + u[j][k-1]);
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

template<typename C, typename T> Chain<Vector<T,2>> Diffusion2D_MonteCarlo(T t, T dt, T d[2], T w[2], Delegate<C,T> step) {
    T val;
    unsigned int n = t / dt;
    T dx = sqrt(2*dt);
    uniform_real_distribution<T> pdf(0.0,1.0);              // Distribution for accepting a move
    default_random_engine rng;                              // Set RNG seeding value
    rng.seed(chrono::high_resolution_clock::now().time_since_epoch().count());

    Chain<Vector<T,2>> particles, *p;
    T dw = w[1]-w[0];
    Vector<T,2> vec = {0,0};
    vec.e[1] = w[0] + pdf(rng)*dw;                          // Initial position
    particles.Add(vec);                                     // of seeding particles

    for(int i = 0; i < n; i++) {                            // Loop of timesteps
        p = &particles;
        for(int j = 0; j < particles.N; j++) {              // Loop of particles
            p = p->next;                                    // Next particle
            val = pdf(rng);                                 // Random number
            if(val <= 0.25) {                               // Sampling rule: x direction

                val = p->e.e[0] - dx * step();              // Calculate backward move
                if(val > 0)                                 // Valid move
                    p->e.e[0] = val;

            } else if(val <= 0.5) {

                // Particle at the source
                if(p->e.e[0] == 0 && p->e.e[1] >= w[0] && p->e.e[1] <= w[1]) {
                    vec.e[1] = w[0] + pdf(rng)*dw;
                    p->prev->Add(vec);                      // Add new particle at prev
                    j++;
                }
                p->e.e[0] += dx * step();                   // Forward move particle
                if(p->e.e[0] >= d[0]) {                     // Remove particle
                    p = p->prev;                            // into the postsynaptic
                    p->Remove();
                    j--;
                }
            } else if(val <= 0.75) {                        // y direction

                if(p->e.e[0]) {
                    p->e.e[1] -= dx * step();               // Calculate backward move
                    if(p->e.e[1] <= 0) {                    // Remove particle
                        p = p->prev;                        // outside cleft
                        p->Remove();
                        j--;
                    }
                }
            } else {

                if(p->e.e[0]) {
                    p->e.e[1] += dx * step();               // Calculate backward move
                    if(p->e.e[1] >= d[1]) {                 // Remove particle
                        p = p->prev;                        // outside cleft
                        p->Remove();
                        j--;
                    }
                }
            }
        }
    }

    return particles;
}

template<typename T> Vector<unsigned int,2> Diffusion2D_Source(T d, T w[2], unsigned int n) {
    Vector<unsigned int,2> m;
    m.e[0] = w[0]/d *(n-1) + (T)0.5;
    m.e[1] = w[1]/d *(n-1) + (T)0.5;
    return m;
}

template<typename T> T StepUniform() {
    return (T)1;
}
template<typename T> T StepGaussian() {

}
