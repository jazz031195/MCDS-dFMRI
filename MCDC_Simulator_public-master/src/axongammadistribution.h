//!  CylinderGammaDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from a Gamma distribution of radiis placed in
*              a single voxel structure.
*   \author    Jasmine Nguyen-Duc
*   \date      Septembre 2022
*   \version   0.2
=================================================================================================*/

#ifndef AXONGAMMADISTRIBUTION_H
#define AXONGAMMADISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "Axon.h"
#include "dynamic_sphere.h"



class AxonGammaDistribution 
{
public:

    double dyn_perc;                                /*!< Percentage of dynamic cylinders that swell                                 */ 
    double volume_inc_perc;                         /*!< Percentage of Volume increase of cylinders that start to swell       */                                       
    std::vector<Axon> axons;                        /*!< Axon vector                                                            */
    unsigned num_obstacles;                         /*!< number of cylnders fit inside the substrate                                */
    double alpha;                                   /*!< alpha coefficient of the Gamma distribution                                */
    double beta;                                    /*!< beta coefficient of the gamma distribution                                 */
    double icvf;                                    /*!< Achieved intra-celular volum fraction in the substrate                     */
    float min_radius;                                /*!< Minimum radius to be sampled from the gamma distribution                  */
    double icvf_current;
    
    Eigen::Vector3d min_limits;                     /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits;                     /*!< voxel max limits (if any)                                                  */
    
    bool active_state;
    bool tortuous;

    double duration; 
    std::vector<double> tortuosities;
    double c2;                                      /*!< ODF                                               */
    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    AxonGammaDistribution (){}

    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    AxonGammaDistribution(double, double, unsigned, double, double,double,Eigen::Vector3d &,Eigen::Vector3d &, float min_radius = 0.001, bool active_state = false, double c2 = 1.0, bool tortuous = false);
     
     /*!
     *  \brief Shows a small histogram of the gamma distribution
    */
    void displayGammaDistribution();
    /*!
     *  \brief Samples and constructs a Gamma distribution
    */
    void createGammaSubstrate(ostream& out);


    /*!
     *  \brief Prints the cylinders positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out);
    std::vector<Dynamic_Sphere> GrowAxon(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out);
    bool check_borders(Eigen::Vector3d pos, double distance_to_border);
    bool isSphereColliding(Dynamic_Sphere sph);
    bool find_next_center(Axon ax, Dynamic_Sphere& s, double dist_, double& rad, Eigen::Vector3d& new_pos, Eigen::Vector3d& prev_pos, int axon_id, ostream& out);
    void fiber_collapse(Eigen::Vector3d& new_pos, Eigen::Vector3d& prev_pos, std::vector<Eigen::Vector3d>& centers, int& fibre_collapsed_nbr, ostream& out);
    bool fill_spheres_in_between(Axon ax, Dynamic_Sphere& added_sphere, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii, ostream& out);
    void shrink_sphere_rad(double& rad, double axon_rad, double& shrink_perc, ostream& out);
    void find_target_point (double c2, double radius, Eigen::Vector3d& initial_point , Eigen::Vector3d& target_point);

private:

    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added cylinders.
     *  \param cylinders List of included cylinders.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
    */
    double  computeICVF(std::vector<Axon> &axons, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);

    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);


};

#endif // AXONGAMMADISTRIBUTION_H
