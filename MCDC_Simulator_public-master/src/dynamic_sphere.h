//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere obstacle class derived from an Obstacle. Defines a analyitical sphere of radius
*   r and centered in center
*   \author    Jonathan Rafael
*   \date      2020
*   \version   1.5
=================================================================================================*/

#ifndef DYN_SPHERE_H
#define DYN_SPHERE_H

#include "sphere.h"

class Dynamic_Sphere : public Sphere
{
public:

    static int count;
    bool swell;
    double volume_inc_perc;
    int ax_id;
    double max_radius;
    double min_radius;
    bool active_state;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Dynamic_Sphere(){id = count++;}
    /*!
     *  \brief Default destructor. Does nothing
     */
    ~Dynamic_Sphere(){count--;}

    /*!
     *  \param center Sphere origin
     *  \param radius Sphere's radius
     *  \param scale  overall scale for when reading files.
     *  \brief Initialize everything.
     */
    Dynamic_Sphere(Eigen::Vector3d center_, double min_radius_, double volume_inc_perc_, bool swell_, int ax_id_,double scale, bool active_state_):
    swell(swell_),
    volume_inc_perc(volume_inc_perc_), 
    ax_id(ax_id_),
    min_radius(min_radius_*scale), 
     active_state(active_state_){
        radius = min_radius;
        if (swell){
            max_radius = sqrt(1+volume_inc_perc)*radius;
        }
        else {
            max_radius = min_radius;
        }
        if (active_state){
            radius = max_radius;
        }
        center = center_*scale;
        id = count++;
    }

    /*!
     *  \brief constrcutor by copy
     */
    Dynamic_Sphere(Dynamic_Sphere const &sph);
    Dynamic_Sphere(Sphere const& s);
    Dynamic_Sphere(Eigen::Vector3d soma_center, double soma_radius);

    /*! \fn  minDistance
     *  \param walker, Walker instance in the simulation.
     *  \brief Returns the minimum distance from the walker to the Sphere. Used to set the reachable
     *  Spheres that a given walker can reach.
     */
    bool isInside(Walker &w);
    bool isInside(Eigen::Vector3d pos, double distance_to_be_inside) const;
    bool distSmallerThan(Eigen::Vector3d pos, double distance);
    void set_center(Eigen::Vector3d center_);

};  

#endif // Sphere_H
