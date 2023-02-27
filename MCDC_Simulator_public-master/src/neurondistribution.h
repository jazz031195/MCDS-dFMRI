//!  NeuronDistribution Class ====================================================================/
/*!
*   \details   Class to construct a substrate of neurons placed in a single voxel structure.
*   \author    Inès de Riedmatten
*   \date      Janvier 2023
*   \version   0.2
=================================================================================================*/

#ifndef NEURONDISTRIBUTION_H
#define NEURONDISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "Neuron.h"
#include "dynamic_sphere.h"
#include <tuple>


class NeuronDistribution 
{
public:
    std::vector<Neuron> neurons;                    /*!< All neurons of the simulation                                              */
    Eigen::Vector3d min_limits_vx;                  /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits_vx;                  /*!< voxel max limits (if any)                                                  */

    NeuronDistribution(){}
    /**
     *  @brief Constructor.
     *  @param num_obstacles              int, number of obstacles in the distribution.
     *  @param icvf                    double, intra compartment volume fraction.
     *  @param min_limits_vx_ Eigen::Vector3d, minimum limit of the simulation voxel.
     *  @param max_limits_vx_ Eigen::Vector3d, maximum limit of the simulation voxel.
    */
    NeuronDistribution(int const& num_obstacles, double const& icvf, Eigen::Vector3d const& min_limits_vx_, Eigen::Vector3d const& max_limits_vx_);
    /**
     *  @brief Populate the simulation voxel with Neurons.
    */ 
    void createSubstrate();
    /**
     *  @brief Prints the neurons positions in a file or output stream.
     *  @param out std::ostream where to write the info.
    */
    void printSubstrate(std::ostream &out) const;
    
private:

    int num_obstacles;                              /*!< number of neurons to fit inside the substrate                              */
    double icvf;                                    /*!< Achieved intra-celular volum fraction in the substrate                     */
    struct projection_pt{                           /*!< Structure to calculate the projection of a sphere                          */
        double position;
        int neuron_id;
        int sphere_id;
    };
    std::vector<projection_pt> projections_x;       /*!< Projections of all dendrites spheres on x axis                             */
    std::vector<projection_pt> projections_y;       /*!< Projections of all dendrites spheres on y axis                             */
    std::vector<projection_pt> projections_z;       /*!< Projections of all dendrites spheres on z axis                             */
    /**
     *  Compute the intracompartment volume fraction (ICVF) of the substrate. 
     *  @return ICVF double.
    */
    double computeICVF() const;
    /**
     *  Compute the max_limits_vx based on the target icvf and the radiis of the axons 
    */
    // void computeMinimalSize(std::vector<double> const& radiis, double &icvf_, Eigen::Vector3d &l) const;
    /**
     * Check if a sphere sph is colliding with this.
     * 
     * @param sph Dynamic_Sphere, sphere to test.
    */
    bool isSphereColliding(Dynamic_Sphere const& sph);
    /**
     * Check if a sphere characterized by its center sphere_center and 
     * its radius sphere_radius is colliding with this.
     * 
     * @param sphere_center Eigen::Vector3d, center of the sphere.
     * @param sphere_radius double         , radius of the sphere.
    */
    bool isSphereColliding(Eigen::Vector3d const& sphere_center, double const& sphere_radius);
    /**
     * Check if a position pos is inside the simulation voxel
     * 
     * @param pos       Eigen::Vector3d, position to check.
     * @param distance_to_border double, minimal distance tolerated from the 
     *                                   borders (e.g. to give space for radius).
     * @return bool, true if the position is inside.
    */
    bool isInVoxel(Eigen::Vector3d const& pos, double const& distance_to_border) const;
    /**
     * Grows the dendrites from a neuron.
     * 
     * @param neuron Neuron.
    */
    void growDendrites(Neuron& neuron);

};

#endif // NEURONDISTRIBUTION_H
