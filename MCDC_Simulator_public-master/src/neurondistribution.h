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


class NeuronDistribution 
{
public:

    std::vector<Neuron> neurons;                    /*!< All neurons of the simulation                                                           */
    unsigned num_obstacles;                         /*!< number of neurons to fit inside the substrate                                */
    double icvf;                                    /*!< Achieved intra-celular volum fraction in the substrate                     */
    Eigen::Vector3d min_limits_vx;                  /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits_vx;                  /*!< voxel max limits (if any)                                                  */
    
    struct projection_pt{
        double position;
        int neuron_id;
        int sphere_id;
    };
    std::vector<projection_pt> projections_x;
    std::vector<projection_pt> projections_y;
    std::vector<projection_pt> projections_z;


    NeuronDistribution (){}

    NeuronDistribution(unsigned num_obstacles, double icvf, Eigen::Vector3d & min_limits, Eigen::Vector3d & max_limits);
     
    void createSubstrate();


    /*!
     *  \brief Prints the cylinders positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out);
    bool isSphereColliding(Dynamic_Sphere sph, double distance_to_be_inside, int axon_id, int sph_id, ostream& out);
    void add_projection(Axon ax, int ax_index, double distance_to_be_inside, ostream& out);
    std::vector<projection_pt> find_collisions(projection_pt proj_on_axis_min, projection_pt proj_on_axis_max,std::vector<projection_pt> projections_on_axis, ostream& out);
    bool isColliding(Axon ax,  double distance_to_be_inside, int axon_id, ostream& out);
    bool search_for_sphere(std::vector<projection_pt> spheres_, projection_pt s);
    bool isInVoxel(Eigen::Vector3d pos, double distance_to_border);
    void growDendrites(Neuron& neuron, int neuron_id);

private:

    /*!
     *  \brief 
     *  \param 
     *  \param 
    */
    double computeICVF(std::vector<Axon> &axons, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);
    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);


};

#endif // NEURONDISTRIBUTION_H
