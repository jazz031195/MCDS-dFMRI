
#ifndef GROWAXONS_H
#define GROWAXONS_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "Axon.h"
#include "dynamic_sphere.h"




class Growth
{
public:

    std::vector<Axon> env_axons; 
    Axon* axon_to_grow;                       /*!< Axon vector                                                            */
    std::vector<Axon> twin_axons;                   /*!< Twins of Axon vector that touches border                                             */
    Eigen::Vector3d small_voxel_size;
    bool tortuous;
    bool success;
    bool finished;

    Growth(){};
    ~Growth();

    Growth (Axon*,  std::vector<Axon>, Eigen::Vector3d, bool);


    void createTwinSpheres(std::vector<Dynamic_Sphere>& twin_spheres, Dynamic_Sphere s, Eigen::Vector2d twin_delta_pos);
    void GrowInParallel(const Eigen::Vector3d destination);
    std::vector<Dynamic_Sphere> GrowAxon(Axon* ax, Eigen::Vector3d destination, std::vector<Eigen::Vector2d>& all_twin_delta_pos, ostream& out);
    void find_shrinking_dichotomy(double& rad, double axon_rad, double min_shrink, double max_shrink, double& half_shrink);
    void add_next_sphere(Dynamic_Sphere added_sphere, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii);
    void fill_with_spheres(Axon* ax, std::vector<Dynamic_Sphere>& spheres_to_add, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii);
    void shrink_sphere_rad(double& rad, double axon_rad, double& shrink_perc);
    bool find_next_center(Axon* ax, Eigen::Vector3d destination, Dynamic_Sphere& s, vector<Eigen::Vector3d> centers, double dist_, double& rad, ostream& out, bool& collides_with_border, Eigen::Vector2d& twin_delta_pos,  std::vector<Eigen::Vector2d> phi_theta_colliding, int max_tries);
    bool isSphereColliding(Dynamic_Sphere sph);
    bool fiber_collapse(std::vector<Eigen::Vector3d>& centers, int fibre_collapsed_nbr, ostream& out);
    bool check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d& twin_delta_pos);
    void createTwinAxons(Axon* ax, Eigen::Vector2d twin_delta_pos);
    
};

#endif 
