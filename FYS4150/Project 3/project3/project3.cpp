#include <fstream>
#include "Gravitation.h"

#define FLOAT double
#define LIST(TYPE) (TYPE*)(const TYPE[])

int main() {
    ofstream file;
    System<FLOAT,2> SolarSystem((FLOAT)4*M_PI*M_PI/333000);

    FLOAT t = 0;
    SolarSystem.Add("Sun", 333000, LIST(FLOAT){0,0}, LIST(FLOAT){0,0});
    SolarSystem.Add("Earth", 1, LIST(FLOAT){0.707,0.707}, LIST(FLOAT){0,0});
    SolarSystem.GoesAround<0,1, CircularDirection::CCW>("Earth", "Sun");
    SolarSystem.Run<DifferentialType::Verlet>(0.5, 10);
    file.open("Sun_and_Earth.dat");
    SolarSystem.Print(file);

    file.close();
}
