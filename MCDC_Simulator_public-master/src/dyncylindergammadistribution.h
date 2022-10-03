//!  CylinderGammaDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from a Gamma distribution of radiis placed in
*              a single voxel structure.
*   \author    Jasmine Nguyen-Duc
*   \date      Septembre 2022
*   \version   0.2
=================================================================================================*/

#ifndef DYNCYLINDERGAMMADISTRIBUTION_H
#define DYNCYLINDERGAMMADISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "cylinder.h"
#include "dynamic_Cylinder.h"
#include "cylindergammadistribution.h"

class DynCylinderGammaDistribution : public CylinderGammaDistribution
{
public:

    double dyn_perc;                                /*!< Percentage of dynamic cylinders that swell                                 */ 
    double activation_time;                         /*!< Time at which the system becomes active and cylinders start to swell       */                                       
    double volume_inc_perc;                         /*!< Percentage of Volume increase of cylinders that start to swell       */                                       
    std::vector<Dynamic_Cylinder> dyn_cylinders;    /*!< Cylinder vector                                                            */


    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    DynCylinderGammaDistribution(){}

    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    DynCylinderGammaDistribution(double dyn_perc, double activation_time,double volume_inc_perc, unsigned num_obstacles,double alpha, double beta,double icvf,Eigen::Vector3d min_limits,Eigen::Vector3d max_limits, float min_radius = 0.001);

    /*!
     *  \brief Samples and constructs a Gamma distribution
    */
    void createGammaSubstrate();

    void updateGammaSubstrate(unsigned t);

    /*!
     *  \brief Prints the cylinders positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out);

private:

    /*!
     *  \brief Checks for collision between inside a voxel (with periodic boundaries)
     *  \param cyl cylinder to check collision with
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param cylinders_list cylinders already added.
     *  \param min_distance that two cylinders can be close to.
    */
    bool checkForCollition(Dynamic_Cylinder cyl, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector<Dynamic_Cylinder>& dyn_cylinders_list, double &min_distance);

    /*!
     *  \brief Auxiliary function to check the BOundary collision
     *  \param cyl cylinder to check collision with.
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param cylinders_list cylinders already added.
    */
    void checkBoundaryConditions(Dynamic_Cylinder cyl, std::vector<Dynamic_Cylinder>& dyn_cylinders_list, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits);

    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added cylinders.
     *  \param cylinders List of included cylinders.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
    */
    double  computeICVF(std::vector<Dynamic_Cylinder> &dyn_cylinders, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);

    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);


};

#endif // DYNCYLINDERGAMMADISTRIBUTION_H
