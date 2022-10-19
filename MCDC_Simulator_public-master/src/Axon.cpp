#include "Axon.h"
#include "dynamic_sphere.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include "simerrno.h"

using namespace Eigen;
using namespace std;

int Axon::count = 0;
Axon::Axon()
{
    id = count++;
}

Axon::~Axon()
{
    count--;
}

Axon::Axon(const Axon &ax)
{

    swell = ax.swell;
    id = count++;
    spheres = ax.spheres; 
    radii = ax.radii;
    centers = ax.centers; 


}



