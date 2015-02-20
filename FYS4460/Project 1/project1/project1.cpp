/*
 *  FYS4460 - Uordnede systemer og perkolasjon - Project 1
 *
 *  Written by:         Eimund Smestad
 *
 *  08.02.2015
 *
 *  c++11 compiler
 */
#ifndef PROJECT_CPP
#define PROJECT_CPP

#include <string>
#include <time.h>
#include "array.h"
#include "boltzmann.h"
#include "delegate.h"
#include "euclidean.h"
#include "gaussian.h"
#include "LJpotential.h"
#include "pbc.h"
#include "run.h"
#include "timestep.h"
#include "unit.h"
#include "verlet.h"
#include "vmd.h"

template<typename M, typename T, typename X, typename V, typename D, typename B> struct MD_Model : public VMD_Vel<X,V> {
    M m;
    T dt;
    D* step;
    B* boundary;
    MD_Model(const M& m, const T& dt=1, string comment = "")
        : VMD_Vel<X,V>(comment + " " + X::unit + ", " + V::unit),
          m(m), dt(dt), step(nullptr), boundary(nullptr) {
    }
    inline operator size_t() {
        return this->data;
    }

    inline MD_Model& operator=(const VMD_Vel<X,V>& other) {
        static_cast<VMD_Vel<X,V>&>(*this) = other;
        return *this;
    }
    inline operator T&() {
        return dt;
    }
    inline operator T() const {
        return dt;
    }
    inline operator M&() {
        return m;
    }
    inline operator M() const {
        return m;
    }
    template<typename U, typename W> void Step(TimeStep<U,W>& data) {
        static_cast<T&>(data) += dt;
        if(step != nullptr)
            (*step)(dt);
    }
    template<typename U, typename W> void Boundary(TimeStep<U,W>& data) {
        if(boundary != nullptr)
            (*boundary)(static_cast<VMD<X,Vector<V,3>>&>(*data));
    }
};

int main() {

    // Floating point precision
    typedef double FLOAT;

    // Units
    typedef angstrom<FLOAT> x;                  // length
    typedef K<FLOAT> T;                         // temperature
    typedef u<FLOAT> m;                         // mass
    typedef angstrom_div_ps<FLOAT> v;           // velocity
    typedef angstrom_div_ps2<FLOAT> a;          // acceleration
    typedef u_mul_angstrom_div_ps2<FLOAT> f;    // force
    typedef u_mul_angstrom2_div_ps2<FLOAT> e;   // energy
    typedef ps<FLOAT> t;                        // time

    // Declarations
    auto mass = m(39.948);
    auto temp = T(100);
    auto epsilon = T(119.8);
    auto dx = x(5.26);
    auto sigma = x(3.405);
    auto dt = t(0.01);
    auto t_end = t(10);
    auto N = ARRAYLIST(size_t, 4, 4, 4);
    string molecule[4];
    molecule[0] = "Ar";
    molecule[1] = "Ar";
    molecule[2] = "Ar";
    molecule[3] = "Ar";
    ofstream time_file;
    time_file.open("time.dat");
    auto t0 = clock();

    // Boltzmann velocity dirstribution
    Boltzmann<FLOAT> boltz(temp, mass);
    Gaussian<v> gauss(boltz.StandardDeviation());
    auto boltz_vel = delegate(gauss, &Gaussian<v>::Distribution_mu_0);

    // Models
    typedef VMD<x,Vector<v,3>> pos;
    typedef VMD_Vel<x,v> vel;
    typedef Vector<a,3> va;
    typedef Array<void,va> acc;
    typedef PBC<pos,vel,acc,void,x,3>::CX pref;
    typedef PBC<pos,vel,acc,void,x,3>::CA aref;
    typedef LJ_Potential<3,x,a,m,f,e,MinImage<x,3>,void,void> Potential;
    typedef Delegate<Potential,void,pref&,aref&> LJ_Solver;
    typedef PBC<pos,vel,acc,LJ_Solver,x,3> Boundary;
    typedef Delegate<Boundary,void,pos&> PeriodicBoundary;
    typedef Delegate<Boundary,void,pos&,vel&,acc&> List;
    typedef VelocityVerlet<2,pos,vel,va,List> VerletModel;
    typedef Delegate<VerletModel,void,t&> VerletSolver;

    // Molecular dynamic model
    MD_Model<m,t,x,v,VerletSolver,PeriodicBoundary> model(mass, dt, "Argon centered cubic lattice");
    // Lennard-Jones potenial
    Potential potential(sigma, mass, 0);
    auto LJ_acc = delegate(potential, &Potential::Acceleration<pref,aref>);
    // Periodic bounary condition
    Boundary pbc(LJ_acc);
    for(size_t i = 0; i < 3; i++) {
        pbc.origin[i] = -0.25*dx;
        pbc.range[i] = dx*(N[i]-0.25);
    }
    auto periodic_boundary = delegate(pbc, &Boundary::Boundary<pos>);
    model.boundary = &periodic_boundary;
    auto list = delegate(pbc, &Boundary::NeighbourList);
    auto dist_vec = delegate(static_cast<MinImage<x,3>&>(pbc), &MinImage<x,3>::Distance<x,x>);
    auto dist = delegate(&Euclidean::Length<Vector<x,3>>);
    auto unit_vec = delegate(&Euclidean::UnitVector<x,x,3>);
    potential.dist_vec = &dist_vec;
    potential.dist = &dist;
    potential.unit_vec = &unit_vec;
    // Verlet differential solver
    VerletModel verlet(static_cast<pos&>(model), static_cast<vel&>(model), list);
    auto verlet_solver = delegate(verlet, &VerletModel::Solve<t>);
    model.step = &verlet_solver;

    // Lattice configuration
    model.LatticeCenteredCubic(molecule, dx, N);
    model.SetVelocityDistribution(boltz_vel);
    auto data = static_cast<vel>(model);
    TimeStep<t,decltype(model)> time(model);

    // Task a and b
    /*time.Open("b.xyz");
    time.Print();
    time.Close();*/

    // Task c
    auto step = delegate(model, &decltype(model)::Step<t,decltype(model)>);
    auto Time = time;
    /*Time.Open("c.xyz");
    Run(Time, t_end, step);
    Time.Close();*/

    // Task d
    auto boundary = delegate(model, &decltype(model)::Boundary<t,decltype(model)>);
    /**time.data = data;
    Time = time;
    Time.Open("d.xyz");
    Run(Time, t_end, step, boundary);
    Time.Close();*/

    // Task g
    potential = epsilon * Boltzmann<FLOAT>::kB;
    /**time.data = data;
    Time = time;
    Time.Open("g.xyz");
    verlet = 0;         // Initialize verlet solver
    t0 = clock();
    Run(Time, t_end, step, boundary);
    time_file << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    Time.Close();*/

    // Task h
    for(size_t i = 0; i < 3; i++)
        pbc[i] = (pbc.range[i]-pbc.origin[i])/(3*sigma);
    *time.data = data;
    Time = time;
    Time.Open("h.xyz");
    verlet = 0;         // Initialize verlet solver
    t0 = clock();
    Run(Time, t_end, step, boundary);
    time_file << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\ ";
    Time.Close();

    time_file.close();

    return 0;
}

#endif // PROJECT_CPP
