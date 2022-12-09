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
#include "simerrno.h"
using namespace std;

/// @brief 
class Axon : public Obstacle
{
public:

    static int count;
    std::vector<Dynamic_Sphere> spheres; 
    double radius;
    double max_radius;
    bool swell;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    double volume_inc_perc;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();


    Axon(double radius_, Eigen::Vector3d begin_,Eigen::Vector3d end_, double volume_inc_perc_, bool swell_ = false, double scale = 1):radius(radius_*scale),begin(begin_*scale),end(end_*scale),volume_inc_perc(volume_inc_perc_), swell(swell_){

        std::vector<Eigen::Vector3d> centers; 
        Eigen::Vector3d squeletton = (end-begin).normalized();
        Eigen::Vector3d pos = begin;
        spheres.clear();


        while ((pos[0] != end[0]) && (pos[1] != end[1])&& (pos[2] != end[2]) ){
            double dist_ = radius/2;
            pos[0] = pos[0] + squeletton[0]*dist_;
            pos[1] = pos[1] + squeletton[1]*dist_;
            pos[2] = pos[2] + squeletton[2]*dist_;
            centers.push_back(pos);
        }


        for (unsigned i=0; i< centers.size(); ++i){
            Dynamic_Sphere sphere(centers[i], radius,volume_inc_perc, swell);  
            spheres.push_back(sphere);
        }
        id = count++;

        max_radius = sqrt(1+volume_inc_perc)*radius;
    }
    Axon(Axon const &ax);

    std::tuple<bool, std::vector<Dynamic_Sphere>> IsInside(Walker &w, double distance_to_be_inside);

    Dynamic_Sphere closestSphere(Walker &w);

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    Dynamic_Sphere closestSphere(Eigen::Vector3d pos);

    std::tuple<bool, std::vector<Dynamic_Sphere>> IsInside(Eigen::Vector3d pos, double distance_to_be_inside);

};

#endif // AXON_H
