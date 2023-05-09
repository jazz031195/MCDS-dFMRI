//!  Cylinder Obstacle Derived Class =============================================================/
/*!
*   \details   Cylinder class derived from an Obstacle. Defines infinite long cylinders
*              in the direction set by P,Q.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2022 
*   \version   1.42
=================================================================================================*/


#ifndef AXON_H
#define AXON_H

#include "dynamic_sphere.h"
#include "swipeprune.h"
#include "simerrno.h"
#include "obstacle.h"
#include <vector>

using namespace std;

/// @brief 
class Axon : public Obstacle
{
public:
    int id;
    std::vector<Dynamic_Sphere> spheres; 
    double radius;
    bool swell;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    double volume_inc_perc;
    Projections projections;
    Projections projections_max;
    double min_radius;


    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon(){};

    ~Axon(){};


    Axon(int id_, double radius_,  Eigen::Vector3d begin_,Eigen::Vector3d end_, double volume_inc_perc_, bool swell_ , int scale = 1){

        id = id_;
        begin = begin_*scale;
        end = end_*scale;
        volume_inc_perc = volume_inc_perc_;
        swell = swell_;

        min_radius = radius = radius_*scale;
        spheres.clear();

        if (swell){ 
            radius = sqrt(1+volume_inc_perc)*radius;
        }

        //projections.clear_projections();
        //projections_max.clear_projections();

    }
    Axon(Axon const &ax);

    std::vector<double> Distances_to_Spheres(Walker &w);

    std::vector<Dynamic_Sphere>  closestSpheres (Walker &w);

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    std::vector<Dynamic_Sphere>  closestSpheres (Eigen::Vector3d pos);

    std::vector<double> Distances_to_Spheres(Eigen::Vector3d pos);
    std::vector<double> Distances_to_Centers(Eigen::Vector3d pos);
    double minDistanceCenter(Eigen::Vector3d pos);
    double minDistance(Walker &w);

    double minDistance(Eigen::Vector3d pos);

    bool intersection_sphere_vector(double &t1, double &t2,Dynamic_Sphere &s, Eigen::Vector3d &step, double &step_length, Eigen::Vector3d &pos, double &c);
    void set_spheres(std::vector<Dynamic_Sphere> spheres_to_add);
    void add_projection();
    bool isNearAxon(Eigen::Vector3d &position,  double distance_to_be_inside);
    bool isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside,  std::vector<int>& sphere_ids);
    
};


#endif // AXON_H
