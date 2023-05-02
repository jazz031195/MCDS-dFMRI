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
    std::vector <int> proximal_branching; /* Contains the parent axon id as well as the other child id at the proximal branching */
    std::vector <int> distal_branching;   /* Contains the child.ren id at the distal branching */


    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();


    Axon(double const& min_radius_,  Eigen::Vector3d const& begin_, Eigen::Vector3d const& end_, double const& volume_inc_perc_, bool const& active_state_ , bool const& swell_ , double const& scale):
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
        proximal_branching.clear();
        distal_branching.clear();

        id = count++;
    }
    
    Axon(double const& min_radius_,  Eigen::Vector3d const& begin_, Eigen::Vector3d const& end_, double const& volume_inc_perc_, 
         bool const& active_state_, bool const& swell_, double const& scale, std::vector <int> proximal_branching_, std::vector <int> distal_branching_,
         int const& id_):
    Axon(min_radius_, begin_, end_, volume_inc_perc_, active_state_ , swell_ , scale)
    {
        proximal_branching = proximal_branching_;
        distal_branching   = distal_branching_;
        id                 = id_;
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
    bool checkCollision(Walker &walker, Eigen::Vector3d const&step, double const&step_lenght, Collision &colision) const;
   
    std::vector<Dynamic_Sphere>  closestSpheres (Eigen::Vector3d const& pos) const;
    std::vector<Dynamic_Sphere>  closestSpheres (Walker const& w) const;
    int closest_sphere_dichotomy(Walker const& walker) const;
    
    std::vector<double> Distances_to_Spheres(Walker const& w) const;
    std::vector<double> Distances_to_Spheres(Eigen::Vector3d const& pos) const;
    std::vector<double> Distances_to_Centers(Eigen::Vector3d const& pos) const;
    double minDistanceCenter(Eigen::Vector3d const& pos) const;
    double minDistance(Walker const& w) const;
    double minDistance(Eigen::Vector3d const& pos) const;

    bool intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere const&s, Eigen::Vector3d const&step, double const&step_length, Eigen::Vector3d const& pos, bool const& isintra, double &c) const;
    void set_spheres(std::vector<Dynamic_Sphere> const& spheres_to_add, int const& axon_id);
    void add_projection(int const& axon_id);
    bool isNearAxon(Eigen::Vector3d const&position,  double const& distance_to_be_inside) const;
    bool isPosInsideAxon(Eigen::Vector3d const&position,  double const& distance_to_be_inside, bool const& swell_, std::vector<int> sphere_ids) const;
    /* Calculate the volume of the current axon and update its tortuosity */
    double volumeAxon() const;
};


#endif // AXON_H
