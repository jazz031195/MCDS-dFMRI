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
    double tortuosity; /* Total length of an axon divided by the distance between begin and end (as if it was a straight cylinder)*/
    Projections projections;
    Projections projections_max;


    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();


    Axon(double min_radius_,  Eigen::Vector3d begin_, Eigen::Vector3d end_, double volume_inc_perc_,bool active_state_ , bool swell_ , double scale):
    min_radius(min_radius_*scale), 
    swell(swell_), 
    begin(begin_*scale),
    end(end_*scale),
    volume_inc_perc(volume_inc_perc_), 
    active_state(active_state_){

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

        projections.clear_projections();
        projections_max.clear_projections();

        id = count++;
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

    bool intersection_sphere_vector(double &t1, double &t2,Dynamic_Sphere &s, Eigen::Vector3d &step, double &step_length, Eigen::Vector3d &pos, bool isintra, double &c);
    void set_spheres(std::vector<Dynamic_Sphere> spheres_to_add, int axon_id);
    void add_projection(int axon_id);
    bool isNearAxon(Eigen::Vector3d &position,  double distance_to_be_inside);
    bool isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, bool swell_, std::vector<int> sphere_ids);
    int closest_sphere_dichotomy(Walker &walker, double &step_lenght);
    /* Calculate the volume of the current axon and update its tortuosity */
    double volumeAxon();
    bool check_negatif(vector<double> list);
};


#endif // AXON_H
