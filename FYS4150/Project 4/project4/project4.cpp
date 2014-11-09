/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  05.11.2014
 *
 *  c11 compiler
 */

#include "HeatEquation.h"
#include "time.h"
#include <fstream>

#define FLOAT double

int main() {

    auto diffusion = HeatEquation<FLOAT,1>(0, ARRAYLIST(FLOAT,0,0), ARRAYLIST(FLOAT,1,1), ARRAYLIST(int,201,10));

    int N[] = {10};
    FLOAT t[] = {1e-1};
    FLOAT a[] = {0.49};
    ofstream file_fe, file_be, file_cn;
    char filename[50];
    for(int i = 0; i < ARRAY_SIZE(N); i++) {
        sprintf(filename,"FE_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_fe.open(filename);
        sprintf(filename,"BE_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_be.open(filename);
        sprintf(filename,"CN_t%g_a%g_N%d.dat",t[i],a[i],N[i]);
        file_cn.open(filename);

        diffusion.upper[0] = t[i];
        diffusion.n[0] = (int)(t[i]*N[i]*N[i]/a[i]);

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
    file_fe.open("time_forward_euler.dat");
    file_be.open("time_backward_euler.dat");
    file_cn.open("time_crank_nicoloson.dat");


    /*diffusion.upper[0] = 1.0;
    for(int i = 0; i < ARRAY_SIZE(nx); i++) {

        if(i) {
            file_fe << endl;
            file_be << endl;
            file_cn << endl;
        }
        file_fe << nx[i];
        file_be << nx[i];
        file_cn << nx[i];
        for(int j = 0; j < i; j++) {
            file_fe << " & ";
            file_be << " & ";
            file_cn << " & ";
        }

        for(int j = i; j < ARRAY_SIZE(nt); j++) {
            diffusion.n[0] = nt[j];
            diffusion.n[1] = nx[i];

            diffusion.theta = 0.0;
            diffusion.Initialize();
            diffusion.u[0][0] = 1;
            clock_t t0 = clock();
            diffusion.Solve();
            file_fe << " & " << (double)(clock()-t0)/CLOCKS_PER_SEC;

            diffusion.theta = 1.0;
            diffusion.Initialize();
            diffusion.u[0][0] = 1;
            t0 = clock();
            diffusion.Solve();
            file_be << " & " << (double)(clock()-t0)/CLOCKS_PER_SEC;

            diffusion.theta = 0.5;
            diffusion.Initialize();
            diffusion.u[0][0] = 1;
            t0 = clock();
            diffusion.Solve();
            file_cn << " & " << (double)(clock()-t0)/CLOCKS_PER_SEC;
        }

        file_fe << " \\\\";
        file_be << " \\\\";
        file_cn << " \\\\";
    }
    file_fe.close();
    file_be.close();
    file_cn.close();*/

    return 0;
}
