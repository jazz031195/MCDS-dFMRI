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

    int id;
    bool swell;
    double volume_inc_perc;
    int ax_id;
    double min_radius;


    /*!
     *  \brief Default constructor. Does nothing
     */
    Dynamic_Sphere(){}
    /*!
     *  \brief Default destructor. Does nothing
     */
    ~Dynamic_Sphere(){}

    /*!
     *  \param center Sphere origin
     *  \param radius Sphere's radius
     *  \param scale  overall scale for when reading files.
     *  \brief Initialize everything.
     */
    Dynamic_Sphere(Eigen::Vector3d center_, double radius_, double volume_inc_perc_, bool swell_, int ax_id_, int id_, double scale=1):
    id(id_), swell(swell_), volume_inc_perc(volume_inc_perc_), ax_id(ax_id_), min_radius(radius_*scale)
    {
        center = center_*scale;
        radius = radius_*scale;
        if (swell){
            radius = sqrt(1+volume_inc_perc)*radius;
        }
    }

    /*!
     *  \brief constrcutor by copy
     */
    Dynamic_Sphere(Dynamic_Sphere const &sph);

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

private:

    /*! \fn  handleCollition
     *  \param walker, Walker instance in the simulation.
     *  \param collision, Collision instance to save all the information.
     *  \param step, step vector where to move.
     *  \brief Returns true if it was any analytical collision to the infinite plane
     */
    inline bool handleCollition(Walker& walker, Collision &colision, Eigen::Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length);

};  

#endif // Sphere_H
