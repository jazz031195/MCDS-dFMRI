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

    int ub = 7; // upper bound for nb_dendrites
    int lb = 1; // lower bound for nb_dendrites

    std::uniform_int_distribution<std::mt19937::result_type> dist_dendrites(lb, ub);
    // Generate int number in [lb, ub]
    nb_dendrites = dist_dendrites(rng);

    ub = 7; // upper bound for span_radius
    lb = 3; // lower bound for span_radius
    std::uniform_int_distribution<std::mt19937::result_type> dist_span_radius(lb, ub);

    // Generate int number in [lb, ub]
    span_radius = dist_span_radius(rng);

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
    growDendrites();
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
}

bool Neuron::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{

}

void Neuron::growDendrites()
{
    std::vector<Eigen::Vector3d> all_start;
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::normal_distribution<double> distribution(0.0, 1.0);

    for(int i = 0; i < nb_dendrites; ++i)
    {

        double x = distribution(generator); 
        double y = distribution(generator);
        double z = distribution(generator);

        while(x==0 & y==0 & z==0)
        {
            x = distribution(generator);
            y = distribution(generator);
            z = distribution(generator);
        }
        double normalization_factor = sqrt(x*x + y*y + z*z);
        x = x/normalization_factor*soma.radius + soma.center[0];
        y = y/normalization_factor*soma.radius + soma.center[1];
        z = z/normalization_factor*soma.radius + soma.center[2];

        Eigen::Vector3d dendrite_start(x, y, z);

        // If the vector is not already contained in all_start, add it. 
        // Otherwise, decrement i and do one more round
        if(i != 0 & std::count(all_start.begin(), all_start.end(), dendrite_start)){ i--; }
        else
        {
            all_start.push_back(dendrite_start);
            Eigen::Vector3d dendrite_direction = dendrite_start - soma.center;
            dendrite_direction.normalize();
            double sphere_radius = 0.5e-3;
            int nb_spheres = span_radius / (4*sphere_radius); //Let's assume that dendrites have a radius of 0.5microns so far
            
            Eigen::Vector3d begin;
            Axon dendrite(sphere_radius, begin, begin, 0, false, false , 1);
            std::vector<Dynamic_Sphere> spheres_to_add;

            for(int j=0; j < nb_spheres; ++j)
            {
                Eigen::Vector3d center = j*dendrite_direction*sphere_radius/4 + dendrite_start;
                Dynamic_Sphere sphere_to_add(center, sphere_radius, 0, false, j, 1, 0);
                spheres_to_add.push_back(sphere_to_add);
            }
            dendrite.set_spheres(spheres_to_add, i);
            dendrites.push_back(dendrite);
        }
    } 
}
