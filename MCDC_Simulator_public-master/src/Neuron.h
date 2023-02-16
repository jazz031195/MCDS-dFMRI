//!  Class to declare several types of neurons =============================================================/
/*!
*   \details   Class derived from an Obstacle. Defines neurons
*   \author    Inès de Riedmatten
*   \date      January 2023 
*   \version   1.42
=================================================================================================*/


#ifndef NEURON_H
#define NEURON_H

#include "Axon.h"
#include "sphere.h"
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

    static int nb_neurons;                  /* Number of neurons in the simulation*/
    int nb_dendrites;                       /* Number of dendrites */
    double span_radius;                     /* Radius [mm] inside which all the dendrites are contained */
    std::vector<Axon> dendrites;            /* Contains all the dendrites of the neuron*/
    Sphere soma;                            /* Soma of the neuron */

    Neuron();

    ~Neuron();

    Neuron(std::vector<Axon> dendrites, Sphere soma);
    Neuron(std::vector<Axon> dendrites, Eigen::Vector3d soma_center, double soma_radius);
    Neuron(Eigen::Vector3d soma_center, double soma_radius);
    
    /* Copy constructor */
    Neuron(Neuron const &neuron);

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);
    /* Vector of all distances between the position pos and the spheres of the neuron */
    std::vector<double> Distances_to_Spheres(Eigen::Vector3d pos);
    std::vector<double> Distances_to_Spheres(Walker &w);

    /* Minimal distance between a position pos and a neuron */
    double minDistance(Eigen::Vector3d pos);
    double minDistance(Walker &w);

    bool isPosInsideNeuron(Eigen::Vector3d &position,  double barrier_thickness, bool swell_);
    /* Check if the position is in the close vicinity (bounding boxes) of an Neuron 
       Returns a tuple {'dendrite', dendrite_id}, {'soma', soma_id} or {'', -1}
    */
    std::tuple<std::string, int> isNearNeuron(Eigen::Vector3d &position,  double barrier_thickness);
    /* Find the closest sphere between a walker and a neuron */
    std::tuple<std::string, int, int> closest_sphere_dichotomy(Walker &walker, double &step_lenght, double barrier_tickness);
    bool intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere &s, Eigen::Vector3d &step, double &step_length, Eigen::Vector3d &pos, bool isintra, double &c);


};


#endif // NEURON_H
