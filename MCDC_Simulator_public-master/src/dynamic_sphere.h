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

#include "obstacle.h"

class Dynamic_Sphere : public Obstacle
{
public:

    int id;
    Eigen::Vector3d center;    /*!< Cilinder Axis reference Points, P should be the "center"      */
    double radius;             /*!< Radius of the Sphere                                          */
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
    Dynamic_Sphere(int id_, int ax_id_, Eigen::Vector3d center_, double volume_inc_perc_, bool swell_, double min_radius_ = 0.0, double max_radius = 0.0): id(id_), center(center_),  volume_inc_perc(volume_inc_perc_), ax_id(ax_id_),swell(swell_){

        if (min_radius_ == 0.0 && max_radius != 0.0){
            min_radius = max_radius /(sqrt(1+volume_inc_perc));
            if (!swell){
                radius = min_radius;
            }
            else{
                radius = max_radius;
            }
        }
        else if (min_radius_ != 0.0 && max_radius == 0.0){
            min_radius = min_radius_;
            if (swell){
                radius = min_radius_*sqrt(1+volume_inc_perc);
            }
            else{
                radius = min_radius_;
            }
        }
        else{
            std::cout << "Please give a valid radius for axon : " << id << endl;
        }
    }

    /*!
     *  \brief constrcutor by copy
     */
    Dynamic_Sphere(Dynamic_Sphere const &sph);

    /*! \fn  checkCollision
     *  \param walker, Walker instance in the simulation.
     *  \param 3d step. Is assumed to be normalized.
     *  \param step_length, length used as the maximum step collision distance.
     *  \param collision, Collision instance to save the collision (if any) details.
     *  \return true only if there was a Collision::hit status. \see Collision.
     *  \brief Basic collision function. Returns the if there was any collision on against the obstacle.
     */
    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    /*! \fn  minDistance
     *  \param walker, Walker instance in the simulation.
     *  \brief Returns the minimum distance from the walker to the Sphere. Used to set the reachable
     *  Spheres that a given walker can reach.
     */
    double minDistance(Walker &w);
    bool isInside(Walker &w);
    bool isInside(Eigen::Vector3d pos, double distance_to_be_inside);
    bool distSmallerThan(Eigen::Vector3d pos, double distance);
    void set_center(Eigen::Vector3d center_);
    double minDistance(Eigen::Vector3d O);
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
