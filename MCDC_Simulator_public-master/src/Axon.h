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


/// @brief 
class Axon : public Obstacle
{
public:

    static int count;
    std::vector<Dynamic_Sphere> spheres; 
    double radius;
    std::vector<Eigen::Vector3d> centers; 
    bool swell;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    double volume_inc_perc;
    double activation_time;
    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();


    Axon(double radius_, std::vector<Eigen::Vector3d> centers_,Eigen::Vector3d begin_,Eigen::Vector3d end_, double volume_inc_perc_,double activation_time_,bool swell_ = false):radius(radius_),centers(centers_),begin(begin_),end(end_),volume_inc_perc(volume_inc_perc_),activation_time(activation_time_), swell(swell_){


        Eigen::Vector3d squeletton = (end-begin).normalized();
        Eigen::Vector3d pos = begin;
        centers.clear();
        while (pos != end){
            double dist_ = radius/2;
            pos = pos + squeletton*dist_;
            centers.push_back(pos);
        }

        for (unsigned i=0; i< centers.size(); ++i){
            Dynamic_Sphere sphere(centers_[i], radius,volume_inc_perc_,activation_time_, swell_);  
            spheres[i] = sphere;
        }
        id = count++;
    }
    Axon(Axon const &ax);

    std::tuple<bool, Dynamic_Sphere> IsInside(Walker &w);

    Dynamic_Sphere closestSphere(Walker &w);

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    Dynamic_Sphere closestSphere(Eigen::Vector3d pos);

    std::tuple<bool, Dynamic_Sphere> IsInside(Eigen::Vector3d pos);

};

#endif // AXON_H
