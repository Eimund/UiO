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
#include "array.h"
#include "boltzmann.h"
#include "delegate.h"
#include "gaussian.h"
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
    typedef angstrom<FLOAT> x;          // length
    typedef K<FLOAT> T;                 // temperature
    typedef u<FLOAT> m;                 // mass
    typedef angstrom_div_ps<FLOAT> v;   // velocity
    typedef m_div_s2<FLOAT> a;          // acceleration
    typedef ps<FLOAT> t;                // time

    // Declarations
    auto mass = m(39.948);
    auto temp = T(100);
    auto dx = x(5.26);
    auto dt = t(0.1);
    auto t_end = t(100);
    auto N = ARRAYLIST(size_t, 8, 8, 8);
    string molecule[4];
    molecule[0] = "Ar";
    molecule[1] = "Ar";
    molecule[2] = "Ar";
    molecule[3] = "Ar";

    // Boltzmann velocity dirstribution
    Boltzmann<FLOAT> boltz(temp, mass);
    Gaussian<v> gauss(boltz.StandardDeviation());
    auto boltz_vel = delegate(gauss, &Gaussian<v>::Distribution_mu_0);

    // Models
    typedef VMD<x,Vector<v,3>> pos;
    typedef VMD_Vel<x,v> vel;
    typedef VelocityVerlet<2,pos,vel,Vector<a,3>,void> BasicVerletModel;
    typedef Delegate<BasicVerletModel,void,t&> BasicVerletSolver;
    typedef Delegate<PBC<x,3>,void,pos&> PeriodicBoundary;

    // Basic molecular dynamic model
    MD_Model<m,t,x,v,BasicVerletSolver,PeriodicBoundary> basic_model(mass, dt, "Argon centered cubic lattice");
    // Differential solver
    BasicVerletModel basic_verlet(static_cast<pos&>(basic_model), static_cast<vel&>(basic_model));
    auto basic_verlet_solver = delegate(basic_verlet, &decltype(basic_verlet)::Solve<t>);
    basic_model.step = &basic_verlet_solver;
    // Periodic bounary condition
    PBC<x,3> pbc;
    for(size_t i = 0; i < 3; i++)
        pbc.range[i] = dx*(N[i]-0.5);
    auto periodic_boundary = delegate(pbc,&decltype(pbc)::Boundary<pos>);
    basic_model.boundary = &periodic_boundary;
    //basic_model.boundary = &pbc;

    // Lattice configuration
    basic_model.LatticeCenteredCubic(molecule, dx, N);
    basic_model.SetVelocityDistribution(boltz_vel);
    auto data = static_cast<vel>(basic_model);
    TimeStep<t,decltype(basic_model)> basic_time(basic_model);

    // Task a and b
    basic_time.Open("b.xyz");
    basic_time.Print();
    basic_time.Close();

    // Task c
    auto basic_step = delegate(basic_model, &decltype(basic_model)::Step<t,decltype(basic_model)>);
    auto BasicTime = basic_time;
    BasicTime.Open("c.xyz");
    Run(BasicTime, t_end, basic_step);
    BasicTime.Close();

    // Task d
    auto basic_boundary= delegate(basic_model,&decltype(basic_model)::Boundary<t,decltype(basic_model)>);
    *basic_time.data = data;
    BasicTime = basic_time;
    BasicTime.Open("d.xyz");
    Run(BasicTime, t_end, basic_step, basic_boundary);
    BasicTime.Close();

    return 0;
}

#endif // PROJECT_CPP
