
/*
 *  FYS4150 - Computational Physics - Project 5
 *
 *  Written by: Eimund Smestad
 *
 *  14.11.2014
 *
 *  c11 compiler
 */

#include <cmath>
#include "time.h"
#include <fstream>

using namespace std;

#define FLOAT double

template<typename T, unsigned int D> struct Vector {
    T element[D];
    ~Vector() {
    }
    operator T*() {
        return element;
    }
    Vector & operator= (const T other[D]) {
        for(int i = 0; i < D; i++)
            element[i] = other[i];
    }
};

template<typename T> T** ArrayAllocate2(unsigned int n[2]);
template<typename T> void ArrayDeallocate2(T** array, unsigned int n);
template<typename T> T* ArrayRange(T lower, T upper, unsigned int n);
template<typename T> void ArrayToFile(ofstream& file, T* array, unsigned int n);
template<typename T> void ArrayToFile2(ofstream& file, T** array, unsigned int n[2]);
template<typename T> Vector<T,2> Diffusion2D_alpha(T alpha, T d[2], unsigned int n[2]);
template<typename T> T Diffusion2D_deltaT(T alpha, T dx2[2]);
template<typename T> Vector<T,2> Diffusion2D_deltaX2(T d[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Exact(T u0, T D, T t, T d[2], T w[2], T* x[2], unsigned int n[2], unsigned int ne[2]);
template<typename T> T** Diffusion2D_Explicit(T alpha, T t, T d[2], T w[2], unsigned int n[2]);
template<typename T> T** Diffusion2D_Initialize(unsigned int m[2], unsigned int n[2]);
template<typename T> Vector<unsigned int,2> Diffusion2D_Source(T d, T w[2], unsigned int n);

int main() {

    ofstream file;

    FLOAT u0 = 1;
    FLOAT D = 1;
    FLOAT t = 10;
    FLOAT* x[2];
    Vector<unsigned int,2> n = {100,100};
    Vector<unsigned int,2> ne = {100,100};
    Vector<FLOAT,2> d = {1,10};
    Vector<FLOAT,2> w = {7,8};
    x[0] = ArrayRange<FLOAT>(0,d[0],n[0]);
    x[1] = ArrayRange<FLOAT>(0,d[1],n[1]);

    file.open("Exact_2D.dat");
    ArrayToFile(file, x[0], n[0]);
    file << endl;
    ArrayToFile(file, x[1], n[1]);
    file << endl;
    auto u = Diffusion2D_Exact<FLOAT>(u0, D, t, d, w, x, n, ne);
    ArrayToFile2(file, u, n);
    ArrayDeallocate2(u, n[0]);
    file.close();

    file.open("Explicit_2D.dat");
    ArrayToFile(file, x[0], n[0]);
    file << endl;
    ArrayToFile(file, x[1], n[1]);
    file << endl;
    u = Diffusion2D_Explicit<FLOAT>(0.1, t, d, w, n);
    ArrayToFile2(file, u, n);
    ArrayDeallocate2(u, n[0]);
    file.close();

    delete [] x[0];
    delete [] x[1];


    return 0;
}

template<typename T> T** ArrayAllocate2(unsigned int n[2]) {
    T** array = new T*[n[0]];
    for(int i = 0; i < n[0]; i++)
        array[i] = new T[n[1]];
    return array;
}

template<typename T> void ArrayDeallocate2(T** array, unsigned int n) {
    for(int i = 0; i < n; i++)
        delete [] array[i];
    delete [] array;
}

template<typename T> T* ArrayRange(T lower, T upper, unsigned int n) {
    T* array = new T[n];
    if(n) {
        T step = (upper-lower)/(n-1);
        array[0] = lower;
        for(int i = 1; i < n; i++)
            array[i] = array[i-1] + step;
    }
    return array;
}

template<typename T> void ArrayToFile(ofstream& file, T* array, unsigned int n) {
    for(unsigned int i = 0; i < n-1; i++)
        file << array[i] << '\t';
    file << array[n-1];
}

template<typename T> void ArrayToFile2(ofstream& file, T** array, unsigned int n[2]) {
    for(int i = 0; i < n[0]; i++) {
        if(i)
            file << endl;
        ArrayToFile(file, array[i], n[1]);
    }
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
    T** u = ArrayAllocate2<FLOAT>(n);

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
    T a1 = a.element[0];
    T a2 = a.element[1];
    ArrayDeallocate2(v,n[0]);

    return u;
}

template<typename T> T** Diffusion2D_Initialize(unsigned int m[2], unsigned int n[2]) {
    T** u = ArrayAllocate2<FLOAT>(n);

    for(int i = 0; i < n[0]; i++) {
        for(int j = 0; j < n[1]; j++)
            u[i][j] = 0;
    }
    for(int i = m[0]; i <= m[1]; i++)
        u[0][i] = 1;

    return u;
}

template<typename T> Vector<unsigned int,2> Diffusion2D_Source(T d, T w[2], unsigned int n) {
    Vector<unsigned int,2> m;
    m.element[0] = w[0]/d *(n-1) + (T)0.5;
    m.element[1] = w[1]/d *(n-1) + (T)0.5;
    return m;
}
