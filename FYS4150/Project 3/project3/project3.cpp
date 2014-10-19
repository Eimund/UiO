#include <fstream>
#include "Gravitation.h"
#include "time.h"

#define FLOAT double
#define LIST(TYPE) (TYPE*)(const TYPE[])

int main() {
    ofstream file, time;
    time.open("time.dat");

    System<FLOAT,2> SolarSystem((FLOAT)4*M_PI*M_PI/333000);

    //FLOAT t = 0;
    SolarSystem.Add("Sun", 333000, LIST(FLOAT){0,0}, LIST(FLOAT){0,0});
    SolarSystem.Add("Earth", 1, LIST(FLOAT){1,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Earth", "Sun");

    time << 10 << " & ";
    clock_t t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 10);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    file.open("Verlet_Sun_and_Earth_10.dat");
    SolarSystem.Print(file);
    file.close();
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 10);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    file.open("RK4_Sun_and_Earth_10.dat");
    SolarSystem.Print(file);
    file.close();
    time << 20 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 20);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    file.open("Verlet_Sun_and_Earth_20.dat");
    SolarSystem.Print(file);
    file.close();
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 20);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    file.open("RK4_Sun_and_Earth_20.dat");
    SolarSystem.Print(file);
    file.close();
    time << 50 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 50);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    file.open("Verlet_Sun_and_Earth_50.dat");
    SolarSystem.Print(file);
    file.close();
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 50);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    file.open("RK4_Sun_and_Earth_50.dat");
    SolarSystem.Print(file);
    file.close();
    time << 100 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 100);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    file.open("Verlet_Sun_and_Earth_100.dat");
    SolarSystem.Print(file);
    file.close();
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 100);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    file.open("RK4_Sun_and_Earth_100.dat");
    SolarSystem.Print(file);
    file.close();
    time << 1000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 1000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 1000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    time << 100000000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 100000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 100000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;

    time.close();
}
