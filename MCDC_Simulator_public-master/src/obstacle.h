//!  Obstacle Base Class ==============================================================================/
/*!
*   \details   Father class to define the base of any other obstacle (wall or substrate)
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
 =====================================================================================================*/

#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "collision.h"
#include "walker.h"
#include "Eigen/Core"
class Obstacle
{
public:

    int id;                         /*!< Unique id of the simulation                                                */
    int count_perc_crossings;       /*!< Auxiliar value to count the number of percolatin crossings in a simulation */
    double percolation;             /*!< Percolation value between 0 and 1.                                         */
    double T2;                      /*!< T2 decay, not used by default                                              */

    /*! \fn  Obstacle
     *  \brief Default constructor. Does nothing.
     */
    Obstacle();

    /*! \param walker, Walker instance in the simulation.
     *  \param 3d step. Is assumed to be normalized.
     *  \param step_lenght, length used as the maximum step collision distance.
     *  \param colision, Collision instance to save the collision (if any) details.
     *  \return true only if there was a Collision::hit status. \see Collision.
     *  \brief Basic collision function. Returns the if there was any collision on against the obstacle.
     */
    bool checkCollision(Walker& walker, Eigen::Vector3d const& step, double const& step_lenght, Collision& colision);

    /*! Modifies step_dir so that it contains the step_dir after the bounce
     *  \param ray_origin Vector3d, origin of the incident ray.
     *  \param normal Vector3d, normal of the obstacle on which to bounce.
     *  \param step_lenght double, length used as the maximum step collision distance.
     *  \param step_dir Vector3d, step direction - incident ray.
     */
    void elasticBounceAgainsPlane(Eigen::Vector3d const& ray_origin, Eigen::Vector3d const& normal, double const& step_length, Eigen::Vector3d &step_dir) const;

    /*!
     *  \param  walker to find the (closest) distance.
     *  \brief  Returns the minimum distance of collision.
     */
    double minDistance(Walker const& walker) const;

};

#endif // OBSTACLE_H
