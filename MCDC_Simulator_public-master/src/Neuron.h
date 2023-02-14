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
using namespace std;

/// @brief 
class Neuron : public Obstacle
{
public:

    static int nb_neurons;                  /* Number of neurons in the simulation*/
    int nb_dendrites;                       /* Number of dendrites */
    int span_radius;                        /* Radius inside which all the dendrites are contained */
    std::vector<Axon> dendrites;            /* Contains all the dendrites of the neuron*/
    Sphere soma;                            /* Soma of the neuron */

    Neuron();

    ~Neuron();

    Neuron(std::vector<Axon> dendrites, Sphere soma);
    Neuron(std::vector<Axon> dendrites, Eigen::Vector3d soma_center, double soma_radius);
    Neuron(Eigen::Vector3d soma_center, double soma_radius);
    
    /* Copy constructor */
    Neuron(Neuron const &neuron);

    void growDendrites();
    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

};


#endif // NEURON_H
