/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  05.11.2014
 *
 *  c11 compiler
 */

#include "Array.h"
#include "HeatEquation.h"
#include <cmath>
#include "time.h"
#include <fstream>

#define FLOAT double

FLOAT* Exact(FLOAT t, int n, int ne);

int main() {

    auto diffusion = HeatEquation<FLOAT,1>(0, ARRAYLIST(FLOAT,0,0), ARRAYLIST(FLOAT,1,1), ARRAYLIST(int,201,10));

    int N[] = {10,10,10,100,100};
    int ne[] = {100,100,100,100,100};
    FLOAT t[] = {1e-2,1,1e-2,1e-2,1e-2};
    FLOAT a[] = {0.49,0.49,0.1,0.49,0.51};
    ofstream file_e, file_fe, file_be, file_cn;
    char filename[50];
    for(int i = 0; i < ARRAY_SIZE(N); i++) {
        sprintf(filename,"Exact_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_e.open(filename);
        sprintf(filename,"FE_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_fe.open(filename);
        sprintf(filename,"BE_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_be.open(filename);
        sprintf(filename,"CN_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_cn.open(filename);

        diffusion.upper[0] = t[i];
        diffusion.n[0] = (int)(t[i]*N[i]*N[i]/a[i]);
        diffusion.n[1] = N[i];

        ArrayToFile(file_e, & &diffusion.X[1]);
        file_e << endl;
        auto u = Exact(t[i],N[i],ne[i]);
        ArrayToFile(file_e, u, N[i]);
        delete [] u;

        diffusion.theta = 0.0;
        diffusion.Initialize();
        diffusion.u[0][0] = 1.0;
        diffusion.Solve();
        diffusion.Print(file_fe);

        diffusion.theta = 1.0;
        diffusion.Initialize();
        diffusion.u[0][0] = 1.0;
        diffusion.Solve();
        diffusion.Print(file_be);

        diffusion.theta = 0.5;
        diffusion.Initialize();
        diffusion.u[0][0] = 1.0;
        diffusion.Solve();
        diffusion.Print(file_cn);

        file_e.close();
        file_fe.close();
        file_be.close();
        file_cn.close();
    }

    int nx[] = {10,100,1000};
    int nt[] = {201,20001,2000001};
    int n[2];
    file_fe.open("time_header.dat");
    file_fe << "$n_x \\backslash n_t$";
    for(int i = 0; i < ARRAY_SIZE(nt); i++)
        file_fe << " & " << nt[i];
    file_fe << " \\\\";
    file_fe.close();
    ofstream file_fe_e, file_be_e, file_cn_e;
    file_fe.open("time_forward_euler.dat");
    file_be.open("time_backward_euler.dat");
    file_cn.open("time_crank_nicoloson.dat");
    file_fe_e.open("time_forward_euler_error.dat");
    file_be_e.open("time_backward_euler_error.dat");
    file_cn_e.open("time_crank_nicoloson_error.dat");


    diffusion.upper[0] = 1e-2;
    for(int i = 0; i < ARRAY_SIZE(nx); i++) {

        if(i) {
            file_fe << endl;
            file_be << endl;
            file_cn << endl;
            file_fe_e << endl;
            file_be_e << endl;
            file_cn_e << endl;
        }
        file_fe << nx[i];
        file_be << nx[i];
        file_cn << nx[i];
        file_fe_e << nx[i];
        file_be_e << nx[i];
        file_cn_e << nx[i];
        /*for(int j = 0; j < i; j++) {
            file_fe << " & ";
            file_be << " & ";
            file_cn << " & ";
            file_fe_e << " & ";
            file_be_e << " & ";
            file_cn_e << " & ";
        }*/

        for(int j = 0; j < ARRAY_SIZE(nt); j++) {
            diffusion.n[0] = nt[j];
            diffusion.n[1] = nx[i];

            auto u = Exact(1e-2,nx[i],100);

            diffusion.theta = 0.0;
            diffusion.Initialize();
            diffusion.u[0][0] = 1;
            clock_t t0 = clock();
            diffusion.Solve();
            file_fe << " & " << (double)(clock()-t0)/CLOCKS_PER_SEC;
            file_fe_e <<  " & " << ArrayRelativeError(u, diffusion.u[0], nx[i]-1);

            diffusion.theta = 1.0;
            diffusion.Initialize();
            diffusion.u[0][0] = 1;
            t0 = clock();
            diffusion.Solve();
            file_be << " & " << (double)(clock()-t0)/CLOCKS_PER_SEC;
            file_be_e <<  " & " << ArrayRelativeError(u, diffusion.u[0], nx[i]-1);

            diffusion.theta = 0.5;
            diffusion.Initialize();
            diffusion.u[0][0] = 1;
            t0 = clock();
            diffusion.Solve();
            file_cn << " & " << (double)(clock()-t0)/CLOCKS_PER_SEC;
            file_cn_e <<  " & " << ArrayRelativeError(u, diffusion.u[0], nx[i]-1);

            delete [] u;
        }

        file_fe << " \\\\";
        file_be << " \\\\";
        file_cn << " \\\\";
        file_fe_e << " \\\\";
        file_be_e << " \\\\";
        file_cn_e << " \\\\";
    }
    file_fe.close();
    file_be.close();
    file_cn.close();
    file_fe_e.close();
    file_be_e.close();
    file_cn_e.close();

    return 0;
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
