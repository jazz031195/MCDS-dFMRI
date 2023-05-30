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
    std::vector <int> proximal_branching; /* Contains the parent axon id as well as the other child id at the proximal branching */
    std::vector <int> distal_branching;   /* Contains the child.ren id at the distal branching */
    double min_radius;


    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon(){};

    ~Axon(){};


    Axon(int const& id_, double const& radius_,  Eigen::Vector3d const& begin_,Eigen::Vector3d const& end_, double const& volume_inc_perc_,
         bool const&swell_, double const& scale=1):
         radius(radius_*scale), swell(swell_), begin(begin_*scale),end(end_*scale), volume_inc_perc(volume_inc_perc_){
 
        min_radius = radius;
        spheres.clear();

        if (swell){ 
            radius = sqrt(1+volume_inc_perc)*radius;
        }

        // projections.clear_projections();
        // projections_max.clear_projections();
        proximal_branching.clear();
        distal_branching.clear();

    }
    
    Axon(int const& id_, double const& radius_,  Eigen::Vector3d const& begin_, Eigen::Vector3d const& end_, double const& volume_inc_perc_, 
         bool const& swell_, double const& scale, std::vector <int> proximal_branching_, std::vector <int> distal_branching_):
    Axon(id_, radius_, begin_, end_, volume_inc_perc_, swell_ , scale)
    {
        proximal_branching = proximal_branching_;
        distal_branching   = distal_branching_;
    }

    Axon(Axon const &ax);
    /**
     * Check the collision within a dendrite and with the soma of the neuron.
     * 
     * @param walker        Walker.
     * @param step Eigen::Vector3d, direction of the next step.
     * @param step_length   double, length of the step.
     * @param colision   Collision, to handle the eventual collision.
     * @param Soma  Dynamic_Sphere, soma of the neuron.
     * @return true if there is a collision, else false.    
     * */
    bool checkCollision(Walker &walker, Eigen::Vector3d const&step, double const&step_lenght, Collision &colision);
   
    std::vector<Dynamic_Sphere>  closestSpheres(Eigen::Vector3d const& pos) const;
    std::vector<Dynamic_Sphere>  closestSpheres(Walker const& w) const;
    
    std::vector<double> Distances_to_Spheres(Walker const& w) const;
    std::vector<double> Distances_to_Spheres(Eigen::Vector3d const& pos) const;
    std::vector<double> Distances_to_Centers(Eigen::Vector3d const& pos) const;
    double minDistanceCenter(Eigen::Vector3d const& pos) const;
    double minDistance(Walker const& w) const;
    double minDistance(Eigen::Vector3d const& pos) const;


    bool intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere const&s, 
                                    Eigen::Vector3d const&step, double const&step_length,
                                    Eigen::Vector3d const&pos, double &c);
    void set_spheres(std::vector<Dynamic_Sphere>& spheres_to_add);
    void add_projection();
    bool isNearAxon(Eigen::Vector3d const&position,  double const& distance_to_be_inside);
    bool isPosInsideAxon(Eigen::Vector3d const&position,  double const& distance_to_be_inside,  std::vector<int>& sphere_ids);
    /* Calculate the volume of the current axon and update its tortuosity */
    double volumeAxon() const;

};


#endif // AXON_H
