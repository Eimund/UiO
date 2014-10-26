#include <fstream>
#include "Gravitation.h"
#include "time.h"

#define FLOAT double
#define LIST(TYPE) (TYPE*)(const TYPE[])

int main() {
    ofstream file, time;
    time.open("time.dat");

    System<FLOAT,2> SolarSystem((FLOAT)4*M_PI*M_PI/333000);

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
    time << 100000000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(1, 100000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(1, 100000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    time.close();

    SolarSystem.Add("Jupiter", 317.8, LIST(FLOAT){5.204267,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Jupiter", "Sun");
    SolarSystem.ZeroMomentum("Sun");

    SolarSystem.Run<DifferentialType::Verlet>(11.8618, 1000);
    file.open("Verlet_Sun_Earth_Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(11.8618, 1000);
    file.open("RK4_Sun_Earth_Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    auto jupiter = SolarSystem.Find("Jupiter");
    jupiter->m *= 10;
    SolarSystem.ZeroMomentum("Sun");
    SolarSystem.Run<DifferentialType::Verlet>(11.8618, 1000);
    file.open("Verlet_Sun_Earth_10Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(11.8618, 1000);
    file.open("RK4_Sun_Earth_10Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    jupiter->m *= 10;
    SolarSystem.ZeroMomentum("Sun");
    SolarSystem.Run<DifferentialType::Verlet>(11.8618, 1000);
    file.open("Verlet_Sun_Earth_100Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(11.8618, 1000);
    file.open("RK4_Sun_Earth_100Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    jupiter->m *= 10;
    SolarSystem.ZeroMomentum("Sun");
    SolarSystem.Run<DifferentialType::Verlet>(11.8618, 1000);
    file.open("Verlet_Sun_Earth_1000Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    SolarSystem.Run<DifferentialType::RK4>(11.8618, 1000);
    file.open("RK4_Sun_Earth_1000Jupiter_1000.dat");
    SolarSystem.Print(file);
    file.close();
    jupiter->m /= 1000;

    SolarSystem.Add("Mercury", 0.055, LIST(FLOAT){0.387098,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Mercury", "Sun");
    SolarSystem.Add("Venus", 0.815, LIST(FLOAT){0.723327,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Venus", "Sun");
    SolarSystem.Add("Mars", 0.107, LIST(FLOAT){1.523679,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Mars", "Sun");
    SolarSystem.Add("Saturn", 95.152, LIST(FLOAT){9.5820172,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Saturn", "Sun");
    SolarSystem.Add("Uranus", 14.536, LIST(FLOAT){19.189253,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Uranus", "Sun");
    SolarSystem.Add("Neptune", 17.147, LIST(FLOAT){30.070900,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Neptune", "Sun");
    SolarSystem.Add("Pluto", 0.00218, LIST(FLOAT){39.264,0}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Pluto", "Sun");
    SolarSystem.ZeroMomentum("Sun");

    time.open("time2.dat");
    time << 40000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(250, 40000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    file.open("Verlet_SolarSystem_10000.dat");
    SolarSystem.Print(file);
    file.close();
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(250, 40000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    file.open("RK4_SolarSystem_10000.dat");
    SolarSystem.Print(file);
    file.close();
    time << 400000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(250, 400000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(250, 400000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    time << 4000000 << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::Verlet>(250, 4000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
    t0 = clock();
    SolarSystem.Run<DifferentialType::RK4>(250, 4000000);
    time << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
    time.close();

    test;


}
