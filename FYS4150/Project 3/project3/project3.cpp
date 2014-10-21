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
    file.open("RK4_Energy_momentum_100.dat");
    SolarSystem.PrintEnergyMomentum(file);
    file.close();
    time << 1000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 1000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 1000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    auto earth = SolarSystem.Find("Earth");
    auto sun = SolarSystem.Find("Sun");
    earth->v[0][0] = sqrt(2*SolarSystem.G*sun->m/earth->x[0][0]);
    earth->v[1][0] = 0;
    SolarSystem.Run<DifferentialType::RK4>(100, 1000);
    file.open("RK4_Sun_and_Earth_escape_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Earth", "Sun");
    /*time << 100000000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 100000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 100000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;*/

    SolarSystem.Add("Jupiter", 317.8, LIST(FLOAT){5.458104,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Jupiter", "Sun");

    SolarSystem.Run<DifferentialType::Verlet>(13, 1000);
    file.open("Verlet_Sun_Earth_Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(13, 1000);
    file.open("RK4_Sun_Earth_Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    auto jupiter = SolarSystem.Find("Jupiter");
    jupiter->m *= 10;
    SolarSystem.Run<DifferentialType::Verlet>(13, 1000);
    file.open("Verlet_Sun_Earth_10Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(13, 1000);
    file.open("RK4_Sun_Earth_10Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    jupiter->m *= 10;
    SolarSystem.Run<DifferentialType::Verlet>(13, 1000);
    file.open("Verlet_Sun_Earth_100Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(13, 1000);
    file.open("RK4_Sun_Earth_100Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    jupiter->m *= 10;
    SolarSystem.Run<DifferentialType::Verlet>(13, 1000);
    file.open("Verlet_Sun_Earth_1000Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(13, 1000);
    file.open("RK4_Sun_Earth_1000Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    jupiter->m /= 1000;

    time.close();
}
