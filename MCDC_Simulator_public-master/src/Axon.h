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
    double min_radius;
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


    Axon(double min_radius_,  Eigen::Vector3d begin_,Eigen::Vector3d end_, double volume_inc_perc_,bool active_state_ , bool swell_ , double scale):min_radius(min_radius_*scale), begin(begin_*scale),end(end_*scale),volume_inc_perc(volume_inc_perc_), active_state(active_state_), swell(swell_){

        radius = min_radius;
        spheres.clear();

        if (swell){ 
            max_radius = sqrt(1+volume_inc_perc)*radius;
        }
        else{
            max_radius = radius;
        } 

        if (swell && active_state){
            radius = max_radius;
        } 

        //double dist_ = this->radius/4;

        //Eigen::Vector3d first_pos = begin;
        //centers.push_back(first_pos);

        //do{

        //    pos[2] = pos[2] + dist_;
        //    centers.push_back(pos);
        //}while (pos[2] < end[2] - dist_);
        // to make axon period wrt z

        //Eigen::Vector3d last_pos = end;
        //centers.push_back(last_pos);

        //double position;
        //for (unsigned i=0; i< centers.size(); ++i){
        //    Dynamic_Sphere sphere(centers[i], radius,volume_inc_perc,  swell, id, scale);  
        //    this->spheres.push_back(sphere);

        //}
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

    double intersection_sphere_vector(Dynamic_Sphere &s, Eigen::Vector3d &step, double& step_length, Eigen::Vector3d &pos, bool isintra);

    bool isInside(Eigen::Vector3d pos, double distance_to_be_inside);

    void set_spheres(std::vector<Dynamic_Sphere> &spheres_to_add);

};


#endif // AXON_H
