#include "Neuron.h"
#include "dynamic_sphere.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include "simerrno.h"
#include <cmath>

using namespace Eigen;
using namespace std;

int Neuron::nb_neurons = 0;
Neuron::Neuron()
{
    id = nb_neurons++;
    
    std::random_device dev;
    std::mt19937 rng(dev());

    int ub = 10; // upper bound for nb_dendrites
    int lb = 5; // lower bound for nb_dendrites

    std::uniform_int_distribution<std::mt19937::result_type> dist_dendrites(lb, ub);
    // Generate int number in [lb, ub]
    nb_dendrites = dist_dendrites(rng);

    ub = 5; // upper bound for span_radius 
    lb = 2; // lower bound for span_radius 
    std::uniform_int_distribution<std::mt19937::result_type> dist_span_radius(lb, ub);

    // Generate int number in [lb, ub], in [mm]
    span_radius = dist_span_radius(rng);
    span_radius /=10;

}

Neuron::~Neuron()
{
    nb_neurons--;
}

Neuron::Neuron(std::vector<Axon> dendrites_, Sphere soma_) : Neuron()
{
    dendrites = dendrites_;          
    soma = soma_;
}

Neuron::Neuron(Eigen::Vector3d soma_center, double soma_radius=5e-3) : Neuron()
{
    soma = Sphere(soma_center, soma_radius);
}

Neuron::Neuron(std::vector<Axon> dendrites_, Eigen::Vector3d soma_center, double soma_radius=5e-3) : Neuron()
{
    dendrites = dendrites_;          
    soma = Sphere(soma_center, soma_radius);
}

Neuron::Neuron(const Neuron &neuron)
{
    id = nb_neurons++;
    dendrites = neuron.dendrites;           
    soma = neuron.soma;
    nb_dendrites = neuron.nb_dendrites;
    span_radius = neuron.span_radius;
}

bool Neuron::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{

}

