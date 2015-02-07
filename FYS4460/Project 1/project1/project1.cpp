#ifndef PROJECT_CPP
#define PROJECT_CPP

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "boltzmann.h"
#include "gaussian.h"
#include "vector.h"
#include "unit.h"


#define FLOAT double

using namespace std;

template<typename X, typename V, size_t D> struct Particle {
    string label;
    Vector<X,D> x;
    Vector<V,D> v;
    inline Array<void,Ref<X>> Get_x(ArrayLength<void,Particle<X,V,3>> other, size_t k) {
        Array<void,Ref<X>>  data(other);
        for(size_t i = 0; i < other; i++)
            data[i] = other[i].x[k];
        return data;
    }
    inline Array<void,Ref<V>> Get_v(ArrayLength<void,Particle<X,V,3>> other, size_t k) {
        Array<void,Ref<V>>  data(other);
        for(size_t i = 0; i < other; i++)
            data[i] = other[i].v[k];
        return data;
    }
};

template<typename U1, typename U2, typename C> Array<void,Particle<U1,U2,3>> LatticeCenteredCubic(string e[4], U1 b,  Delegate<C,U2> vel, size_t N[3]);
template<typename U1, typename U2> void LatticeToVMD(ofstream& file, ArrayLength<void,Particle<U1,U2,3>> data, string comment);

int main() {
    auto N = ARRAYLIST(size_t, 8, 8, 8);
    string l[4];
    l[0] = "Ar";
    l[1] = "Ar";
    l[2] = "Ar";
    l[3] = "Ar";

    ofstream file;
    file.open("a/Ar.xyz");

    Boltzmann<FLOAT> boltz(K<FLOAT>(100),u<FLOAT>(39.948));
    Gaussian<m_div_s<FLOAT>> vel(boltz.StandardDeviation());
    Delegate<Gaussian<m_div_s<FLOAT>>,m_div_s<FLOAT>> vel_distr(vel, &Gaussian<m_div_s<FLOAT>>::Distribution_mu_0);
    auto Ar = LatticeCenteredCubic(l, angstrom<FLOAT>(5.26), vel_distr, N);
    LatticeToVMD(file, Ar, "Argon centered cubic lattice");
    auto t = Boltzmann<FLOAT>().kB;

    file.close();
    return 0;
}

template<typename U1, typename U2, typename C> Array<void,Particle<U1,U2,3>> LatticeCenteredCubic(string e[4], U1 b,  Delegate<C,U2> vel, size_t N[3]) {
    Array<void,Particle<U1,U2,3>> data(4*N[0]*N[1]*N[2]);
    for(size_t i = 0, j = 0; j < N[0]; j++) {
        for(size_t k = 0; k < N[1]; k++) {
            for(size_t l = 0; l < N[2]; l++) {
                data[i].label = e[0];             // Particle positions
                data[i].x[0] = b*l;
                data[i].x[1] = b*k;
                data[i].x[2] = b*j;
                i++;
                data[i].label = e[1];
                data[i].x[0] = b*(l+0.5);
                data[i].x[1] = b*(k+0.5);
                data[i].x[2] = b*j;
                i++;
                data[i].label = e[2];
                data[i].x[0] = b*l;
                data[i].x[1] = b*(k+0.5);
                data[i].x[2] = b*(j+0.5);
                i++;
                data[i].label = e[3];
                data[i].x[0] = b*(l+0.5);
                data[i].x[1] = b*k;
                data[i].x[2] = b*(j+0.5);
                i++;
            }
        }
    }
    MathArray<void,Ref<U2>> v;
    for(size_t i = 0; i < 3; i++) {         // Loop over vector components
        v = data->Get_v(data,i);
        v.Get(vel);                         // Random generated velocity
        v -= v.Sum()/v;                     // Remove linear momentum
    }
    return data;
}

template<typename U1, typename U2> void LatticeToVMD(ofstream& file, ArrayLength<void,Particle<U1,U2,3>> data, string comment) {
    file << data << endl << comment;
    for(unsigned int i = 0; i < data; i++) {
        file << endl << data[i].label << " ";
        file << data[i].x[0] << " " << data[i].x[1] << " " << data[i].x[2] << " ";
        file << data[i].v[0] << " " << data[i].v[1] << " " << data[i].v[2];
    }
}

#endif // PROJECT_CPP
