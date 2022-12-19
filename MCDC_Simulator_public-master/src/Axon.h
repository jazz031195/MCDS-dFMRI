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
    bool active_state;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();


    Axon(double radius_, Eigen::Vector3d begin_,Eigen::Vector3d end_, double volume_inc_perc_,bool active_state_ = false, bool swell_ = false, double scale = 1):radius(radius_*scale),begin(begin_*scale),end(end_*scale),volume_inc_perc(volume_inc_perc_), active_state(active_state_), swell(swell_){

        std::vector<Eigen::Vector3d> centers; 
        Eigen::Vector3d pos = begin;
        spheres.clear();
        centers.clear();
        double length = (end-begin).norm();

        double dist_ = this->radius/2;

        centers.push_back(pos);
        do{
            //pos[0] = pos[0] + squeletton[0]*dist_;
            //pos[1] = pos[1] + squeletton[1]*dist_;
            pos[2] = pos[2] + dist_;
            centers.push_back(pos);
        }while (pos[2] < length/2);
        // to make axon period wrt z

        Eigen::Vector3d last_pos = end;
        centers.push_back(last_pos);
        do{
            //last_pos[0] = last_pos[0] - squeletton[0]*dist_;
            //last_pos[1] = last_pos[1] - squeletton[1]*dist_;
            last_pos[2] = last_pos[2] - dist_;
            centers.push_back(last_pos);
        }while (last_pos[2] >= length/2);

        max_radius = sqrt(1+volume_inc_perc)*radius;

        if (swell && active_state){
            radius = max_radius;
        } 


        for (unsigned i=0; i< centers.size(); ++i){
            Dynamic_Sphere sphere(centers[i], radius,volume_inc_perc,  swell, id);  
            this->spheres.push_back(sphere);

        }

        id = count++;



    }
    Axon(Axon const &ax);

    std::vector<double> Distances_to_Spheres(Walker &w);

    std::vector<Dynamic_Sphere>  closestSpheres (Walker &w);

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    std::vector<Dynamic_Sphere>  closestSpheres (Eigen::Vector3d pos);

    std::vector<double> Distances_to_Spheres(Eigen::Vector3d pos);

    double minDistance(Walker &w);

    double minDistance(Eigen::Vector3d pos);

    double intersection_sphere_vector(Dynamic_Sphere &s, Eigen::Vector3d &step, Eigen::Vector3d &pos, bool isintra);

    bool isInside(Eigen::Vector3d pos, double distance_to_be_inside);

};


#endif // AXON_H
