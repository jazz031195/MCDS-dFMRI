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
    radius = ax.radius;
    centers = ax.centers; 
    begin = ax.begin;
    end = ax.end;

}

std::tuple<bool, Dynamic_Sphere> Axon::IsInside(Eigen::Vector3d pos){

    for (unsigned i=0; i< spheres.size(); ++i){
        Vector3d m = pos - spheres[i].center;
        double distance_to_sphere = m.norm() - spheres[i].radius;
        if (distance_to_sphere <= 0){
            return std::make_tuple(true, spheres[i]);
            break;
        }
        
    }
    Dynamic_Sphere s;
    return std::make_tuple(true, s);
}

std::tuple<bool, Dynamic_Sphere> Axon::IsInside(Walker &w){

    Vector3d O;
    w.getVoxelPosition(O);
    for (unsigned i=0; i< spheres.size(); ++i){
        Vector3d m = O - spheres[i].center;
        double distance_to_sphere = m.norm() - spheres[i].radius;
        if (distance_to_sphere <= 0){
            return std::make_tuple(true, spheres[i]);
            break;
        }
        
    }
    Dynamic_Sphere s;
    return std::make_tuple(true, s);
}


Dynamic_Sphere Axon::closestSphere(Eigen::Vector3d pos){
    Vector3d O = pos;
    double d = 10;
    Dynamic_Sphere closest_sphere;
    for (unsigned i=0; i< spheres.size(); ++i){
        Vector3d m = O - spheres[i].center;
        double distance_to_sphere = m.norm();
        if (distance_to_sphere <= d){
            closest_sphere = spheres[i];
            d = distance_to_sphere;
        }
    }
    return closest_sphere;
}

Dynamic_Sphere Axon::closestSphere(Walker &w){
    Vector3d O;
    w.getVoxelPosition(O);
    double d = 10;
    Dynamic_Sphere closest_sphere;
    for (unsigned i=0; i< spheres.size(); ++i){
        Vector3d m = O - spheres[i].center;
        double distance_to_sphere = m.norm();
        if (distance_to_sphere <= d){
            closest_sphere = spheres[i];
            d = distance_to_sphere;
        }
    }
    return closest_sphere;
}

bool Axon::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    bool isinside;
    Dynamic_Sphere sphere;
    tie(isinside, sphere) =IsInside(walker);
    if (isinside){
        Vector3d O;
        walker.getVoxelPosition(O);
        Eigen::Vector3d next_pos = O + step;
        bool isinside_;
        Dynamic_Sphere sphere_;
        tie(isinside_, sphere_) =IsInside(next_pos);
        if (isinside_){
            return false;
        }
        else{
            Dynamic_Sphere outer_sphere = sphere;
            return outer_sphere.checkCollision(walker, step, step_lenght, colision);
        }

    }
    else {
        Dynamic_Sphere closest_sphere = closestSphere(walker);
        if (closest_sphere.checkCollision(walker, step, step_lenght, colision)){
            return true;
        }
        else {
            for (unsigned i=0; i< spheres.size(); ++i){
                if (spheres[i].checkCollision(walker, step, step_lenght, colision)){
                    return true;
                    break;
                }
            }
        }
    }
    return false;

}


