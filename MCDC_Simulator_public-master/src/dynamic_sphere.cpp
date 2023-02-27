#include "dynamic_sphere.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "simerrno.h"
using namespace Eigen;
using namespace std;

int Dynamic_Sphere::count = 0;

Dynamic_Sphere::Dynamic_Sphere(Vector3d soma_center, double soma_radius): center(soma_center), radius(soma_radius)
{
    swell        = false;
    active_state = false;
    id     = count++;
}

Dynamic_Sphere::Dynamic_Sphere(Sphere const& s): center(s.center), radius(s.radius)
{                                    
    swell        = false;
    active_state = false;
    id           = count++;
}

Dynamic_Sphere::Dynamic_Sphere(const Dynamic_Sphere &sph)
{
    center = sph.center;
    radius = sph.radius;
    swell  = sph.swell;
    volume_inc_perc = sph.volume_inc_perc; 
    ax_id = sph.ax_id;
    id    = count++;
    active_state = sph.active_state;
    min_radius   = sph.min_radius;
    max_radius   = sph.max_radius;

}

void Dynamic_Sphere::set_center(Eigen::Vector3d center_)
{
    this->center = center_; 

}


bool Dynamic_Sphere::isInside(Walker &w){

    //Minimum distance to the sphere wall.
    double d_ = minDistance(w);
   // return d_>0.0?d_:0.0;
    return d_ <= 0;
}

bool Dynamic_Sphere::isInside(Eigen::Vector3d pos, double distance_to_be_inside){
    double d_ = (pos - this->center).norm();
 
    d_ = d_-this->radius;
    
   // return d_>0.0?d_:0.0;
    return d_ <= distance_to_be_inside;
}

bool Dynamic_Sphere::distSmallerThan(Eigen::Vector3d pos, double distance){
    double d_ = (pos - this->center).norm();
    
    return d_ <= distance;
}
