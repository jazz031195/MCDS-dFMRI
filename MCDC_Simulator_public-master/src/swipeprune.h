//!  Projections Class =============================================================/
/*!
*   \details   Class to construct projections of all objects on each of the axes. This speeds up the collision check 
*   \author    Jasmine Nguyen-Duc
*   \date      January 2023
*   \version   0.2
=================================================================================================*/

#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include "constants.h"
#include "simerrno.h"
#include <iostream>



class Projections
{
public:

    struct projection_pt{
        double position;
        int axon_id;
        int sph_id;
    };
    // each axon is in a rectangle, rectangle defined by a vector of 3 Vector2d : [(Xsmall, Xlarge), (Ysmall, Ylarge), (Zsmall, Zlarge)]
    std::vector<Eigen::Vector2d> axon_projections;
    // projection of spheres in x, y, and z axes
    std::vector<projection_pt> sph_projections_x;
    std::vector<projection_pt> sph_projections_y;
    std::vector<projection_pt> sph_projections_z;
    
    Projections();
    void clear_projections();
    std::vector<projection_pt> find_collisions(projection_pt const& proj_on_axis_min, projection_pt const& proj_on_axis_max,std::vector<projection_pt> const& projections_on_axis, double const& distance_to_be_inside) const;
    std::vector<std::vector<projection_pt>> find_collisions_all_axes(Eigen::Vector3d const& position, double const& rad, int const& ax_id, double const& distance_to_be_inside) const;
    bool isProjInside(std::vector<projection_pt> const& projs, projection_pt const& p) const;
    void append_right_place(projection_pt const& p1, projection_pt const& p2, int const& axis);
};

#endif // PROJECTIONS_H
