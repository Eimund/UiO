
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

template<typename T> T** ArrayAllocate2(int n[2]);
template<typename T> void ArrayDeallocate2(T** array, int n);
template<typename T> T* ArrayRange(T lower, T upper, int n);
template<typename T> void ArrayToFile(ofstream& file, T* array, int n);
template<typename T> void ArrayToFile2(ofstream& file, T** array, int n[2]);
FLOAT* Exact(FLOAT t, int n, int ne);
template<typename T> T** Exact2D(T u0, T D, T t, T d[2], T w[2], T* x[2], int n[2], int ne[2]);

int main() {

    ofstream file;

    FLOAT u0 = 1;
    FLOAT D = 1;
    FLOAT t = 0.01;
    FLOAT* x[2];
    int n[2] = {100,100};
    int ne[2] = {100,100};
    FLOAT d[2] = {1,10};
    FLOAT w[2] = {7,8};
    x[0] = ArrayRange<FLOAT>(0,d[0],n[0]);
    x[1] = ArrayRange<FLOAT>(0,d[1],n[1]);

    file.open("Exact_2D.dat");
    ArrayToFile(file, x[0], n[0]);
    file << endl;
    ArrayToFile(file, x[1], n[1]);
    file << endl;

    auto u = Exact2D(u0, D, t, d, w, x, n, ne);
    ArrayToFile2(file, u, n);
    ArrayDeallocate2(u, n[0]);

    delete [] x[0];
    delete [] x[1];
    file.close();

    return 0;
}

template<typename T> T** ArrayAllocate2(int n[2]) {
    T** array = new T*[n[0]];
    for(int i = 0; i < n[0]; i++)
        array[i] = new T[n[1]];
    return array;
}
template<typename T> void ArrayDeallocate2(T** array, int n) {
    for(int i = 0; i < n; i++)
        delete [] array[i];
    delete [] array;
}
template<typename T> T* ArrayRange(T lower, T upper, int n) {
    T* array = new T[n];
    if(n) {
        T step = (upper-lower)/(n-1);
        array[0] = lower;
        for(int i = 1; i < n; i++)
            array[i] = array[i-1] + step;
    }
    return array;
}
template<typename T> void ArrayToFile(ofstream& file, T* array, int n) {
    for(unsigned int i = 0; i < n-1; i++)
        file << array[i] << '\t';
    file << array[n-1];
}
template<typename T> void ArrayToFile2(ofstream& file, T** array, int n[2]) {
    for(int i = 0; i < n[0]; i++) {
        if(i)
            file << endl;
        ArrayToFile(file, array[i], n[1]);
    }
}

FLOAT* Exact(FLOAT t, int n, int ne) {
    FLOAT* u = new FLOAT[n];
    for(int i = 0; i < n; i++) {
        u[i] = 1.0-(FLOAT)i/(n-1);
        for(int j = 1; j < ne; j++)
            u[i] -= 2.0/(j*M_PI)*sin(M_PI*j*i/(n-1))*exp(-j*j*M_PI*M_PI*t);
    }
    return u;
}
template<typename T> T** Exact2D(T u0, T D, T t, T d[2], T w[2], T* x[2], int n[2], int ne[2]) {
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
