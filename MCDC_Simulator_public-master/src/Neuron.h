//!  Class to declare several types of neurons =============================================================/
/*!
*   \details   Class derived from an Obstacle. Defines neurons
*   \author    Inès de Riedmatten
*   \date      February 2023 
*   \version   1.42
=================================================================================================*/


#ifndef NEURON_H
#define NEURON_H

#include "Axon.h"
#include "dynamic_sphere.h"
#include "simerrno.h"
#include <iostream>
#include <random>
#include <tuple>
#include <string>
using namespace std;

/// @brief 
class Neuron : public Obstacle
{
public:

uint8_t nb_dendrites;                   /* Number of dendrites */
double span_radius;                     /* Radius [mm] inside which all the dendrites are contained */
std::vector<Axon> dendrites;            /* Contains all the dendrites of the neuron*/
Dynamic_Sphere soma;                    /* Soma of the neuron */

/*! Default constructor.*/
Neuron();
/*! Constructor
*  \param dendrites Vector3d<Axon>, dendrites of This.
*  \param soma Dynamic_Sphere     , soma of This.
*/
Neuron(std::vector<Axon> const& dendrites, Dynamic_Sphere const& soma);
/*! Constructor
*  \param dendrites Vector3d<Axon>  , dendrites of This.
*  \param soma_center Vector3d<Axon>, center of the soma.
*  \param soma_radius double        , radius of the soma.
*/
Neuron(std::vector<Axon> const& dendrites, Eigen::Vector3d const& soma_center, double const& soma_radius);
/*! Constructor
*  \param soma_center Vector3d<Axon>, center of the soma.
*  \param soma_radius double        , radius of the soma.
*/
Neuron(Eigen::Vector3d const&soma_center, double const& soma_radius);
/*! Copy constructor 
*  \param neuron Neuron
*/
Neuron(Neuron const& neuron);
/*! Default destructor.*/
~Neuron();
/**
 * Check if there is a collision with walker and this.
 *
 * @param walker             Walker.
 * @param step_dir  Eigen::Vector3d, direction of the step.
 * @param step_lenght        double, length of the step. 
 * @param collision       Collision, class to handle collision. 
 * @return                     bool, true if collision with this.
 */
bool checkCollision(Walker &walker, Eigen::Vector3d const&step_dir, double const&step_lenght, Collision &collision);
/**
 * Calculate if position is inside this.
 *
 * @param position Eigen::Vector3d, position of the walker.
 * @param barrier_thickness double, thickness of the cellular barrier.
 * @param swell_              bool, if this swells or not. 
 * @return                    bool, true if position is inside this.
 */
bool isPosInsideNeuron(Eigen::Vector3d const& position,  double const& barrier_thickness, bool const& swell_);
/**
 * Minimal distance between a position pos and this.
 *
 * @param pos Eigen::Vector3d, position of the walker.
 * @return minimal distance bewteen `pos` and this.
 */
double minDistance(Eigen::Vector3d const& pos) const;
/**
 * Minimal distance between a walker and this.
 *
 * @param walker .
 * @return minimal distance bewteen `walker` and this.
 */
double minDistance(Walker const& walker) const;
/**
 * Add dendrite to dendrites. 
 *
 * @param dendrite_to_add Axon.
 */
void add_dendrite(Axon const& dendrite_to_add);

private:

    static int nb_neurons;                  /* Number of neurons in the simulation*/
    /**
     * Vector of all distances between the position pos and the spheres of the this.
     * 
     * @param pos Eigen::Vector3d, position of the walker.
     * @return std::vector<double>, all the distances. 
     */    std::vector<double> Distances_to_Spheres(Eigen::Vector3d const& pos) const;
    /**
     * Vector of all distances between the walker w and the spheres of the this.
     * 
     * @param w Walker.
     * @return std::vector<double>, all the distances. 
     */
    std::vector<double> Distances_to_Spheres(Walker const& w) const;
    /**
     * Check if the position is in the close vicinity (bounding boxes) of this.
     * @return std::tuple<std::string, int>, {'neuron_part', part_id}
     *                                       neuron_part : "soma", "dendrite" or "none",
     *                                       part_id     : soma_id, dendrite_id or -1.
    */
    std::tuple<std::string, int> isNearNeuron(Eigen::Vector3d const& position,  double const& barrier_thickness) const;
    /**
     *  Find the closest sphere between a walker at walker_pos and the dendrite dendrite_id,
     *  with dichotomy method.
     * 
     * @param dendrite_id            int, index of the dendrite.           
     * @param walker_pos Eigen::Vector3d, position of the walker.
     * @return i                     int, index of the closest sphere.
     */
    int closest_sphere_dichotomy(int const& dendrite_id, Eigen::Vector3d const& walker_pos) const;
    /**
     *  Find the closest sphere between a walker and this.
     * 
     * @param walker            Walker, object from which to calculate the closest sphere 
     *                                  from this.
     * @param barrier_thickness double, thickness of the cellular barrier.
     * @return std::tuple<std::string, int, int>, {neuron_part, part_id, sphere_id},
     *                                            neuron_part : "soma", "dendrite" or "none",
     *                                            part_id     : soma_id, dendrite_id or -1,
     *                                            sphere_id   : 0      , sphere_id   or -1.
     */
    std::tuple<std::string, int, int> closest_sphere(Walker const& walker, double const& barrier_tickness) const;
    /**
     * Calculates if there is/are intersection(s) between the sphere s and a walker
     * starting at traj_orig, with a direction step_dir. 
     * There can be none, one or two intersections.
     * 
     * Taken from : https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
     *
     * @param intercept1         double, distance from traj_orig to the first potential intercept.
     * @param intercept2         double, distance from traj_orig to the second potential intercept.
     * @param s          Dynamic_sphere, sphere with which to calculate the intersection.
     * @param step_length        double, step length of the walker
     * @param traj_orig Eigen::Vector3d, trajectory origin of the walker
     * @param c                  double, ||traj.orig - s.center||² - s.radius²
     */
    bool intersection_sphere_vector(double &intercept1, double &intercept2, Dynamic_Sphere const&s, Eigen::Vector3d const&step_dir, 
                                    double const&step_length, Eigen::Vector3d const&traj_orig, double &c) const;
    /**
     * Assign axons to dendrites. 
     *
     * @param axons std::vector<Axon>, set of axons.
     */
    void set_axons(std::vector<Axon> const& axons);
    /**
     * Generate a random span radius in [mm], in the interval [lower_bound, upper_bound].
     *
     * @param lower_bound int, lower bound of the interval, in [10*mm].
     * @param upper_bound int, upper bound of the interval, in [10*mm].
     */
    void generateSpanRadius(int const& lower_bound=2, int const& upper_bound=5);
};


#endif // NEURON_H
