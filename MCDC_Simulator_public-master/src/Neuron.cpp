#include "Neuron.h"
#include "dynamic_sphere.h"
#include "sphere.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include "simerrno.h"
#include <cmath>
#include <typeinfo>

using namespace Eigen;
using namespace std;

int Neuron::nb_neurons = 0;
Neuron::Neuron()
{
    id = nb_neurons++;
    
    random_device dev;
    mt19937 rng(dev());

    int ub = 10; // upper bound for nb_dendrites
    int lb = 5; // lower bound for nb_dendrites

    uniform_int_distribution<mt19937::result_type> dist_dendrites(lb, ub);
    // Generate int number in [lb, ub]
    nb_dendrites = dist_dendrites(rng);

    // Create a random span radius and set its value to this
    generateSpanRadius();

}

Neuron::~Neuron()
{
    nb_neurons--;
}

Neuron::Neuron(vector<Axon> dendrites_, Dynamic_Sphere soma_) : Neuron()
{
    dendrites    = dendrites_;     
    nb_dendrites = dendrites_.size();
    soma         = soma_;
}

Neuron::Neuron(Vector3d soma_center, double soma_radius=5e-3) : Neuron()
{
    soma = Dynamic_Sphere(soma_center, soma_radius);
}

Neuron::Neuron(vector<Axon> dendrites_, Vector3d soma_center, double soma_radius=5e-3) : Neuron()
{
    dendrites    = dendrites_; 
    nb_dendrites = dendrites_.size();        
    soma         = Dynamic_Sphere(soma_center, soma_radius);
}

Neuron::Neuron(const Neuron &neuron)
{
    id           = nb_neurons++;
    dendrites    = neuron.dendrites;           
    soma         = neuron.soma;
    nb_dendrites = neuron.nb_dendrites;
    span_radius  = neuron.span_radius;
}

double Neuron::minDistance(Walker &w)
{
    vector<double> distances;
    distances.clear();
    distances  = Distances_to_Spheres(w);
    double min = *min_element(begin(distances), end(distances ));

    return min;
}

vector<double> Neuron::Distances_to_Spheres(Walker &w)
{

    Vector3d O;
    w.getVoxelPosition(O);
    return Distances_to_Spheres(O);

}

double Neuron::minDistance(Eigen::Vector3d pos)
{

    vector<double> distances = Distances_to_Spheres(pos);
    double min = *min_element(begin(distances), end(distances ));
    return min;

}

vector<double> Neuron::Distances_to_Spheres(Vector3d pos)
{
    
    vector<double> distances;
    distances.clear();
    
    // First check distance to soma
    Vector3d m = pos - soma.center;
    double distance_to_sphere = m.norm() - soma.radius;
    distances.push_back(distance_to_sphere);

    // Then iterate through dendrites
    for (unsigned i=0; i < nb_dendrites; ++i)
    {
        for (unsigned j=0; j < dendrites[i].spheres.size(); ++j)
        {       
            Vector3d m = pos - dendrites[i].spheres[j].center;
            double distance_to_sphere = m.norm() - dendrites[i].spheres[j].radius;
            distances.push_back(distance_to_sphere);
        }
    }
    return distances;
}

bool Neuron::isPosInsideNeuron(Eigen::Vector3d &position,  double distance_to_be_inside, bool swell_){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other neurons -> check with max_radius so there is room for swelling
    vector<vector<Projections::projection_pt>> coliding_projs;
    bool colliding_all_axes;
    Dynamic_Sphere sphere_ ;
    double rad;
    double rad_;
    // if position is in box with axon inside
    string neuron_part; // "soma", "dendrite" or "none"
    int part_id; // id of the soma or dendrite. -1 if not in neuron
    tie(neuron_part, part_id) = isNearNeuron(position, distance_to_be_inside);
    if(!(neuron_part=="none")){
        if (neuron_part=="dendrite")
        {
            if (dendrites[part_id].isPosInsideAxon(position, distance_to_be_inside, false))
            {
                return true;
            } 
         } 
        else if (neuron_part == "soma")
        {
            // distance from position to soma boundary
            double dis = soma.minDistance(position);
            if( dis <= distance_to_be_inside ){return true;}
        }  
    }
    return false;
} 


tuple<string, int> Neuron::isNearNeuron(Vector3d &position,  double distance_to_be_inside)
{
    // Check soma box
    int count_isnear = 0;
    for (unsigned int axis=0 ; axis < 3 ; ++axis)
    {
        if ((position[axis] >= soma.center[axis] - soma.radius - distance_to_be_inside) && 
            (position[axis] <= soma.center[axis] + soma.radius + distance_to_be_inside))
        {
            ++count_isnear; 
        }
    }
    if (count_isnear == 3){return tuple<string, int>{"soma", soma.id};}
    
    // Check each dendrite's box
    for (unsigned int i=0 ; i < nb_dendrites ; ++i)
    {
        count_isnear = 0;
        for (unsigned int axis=0 ; axis < 3 ; ++axis)
        {
            Vector2d axis_limits = dendrites[i].projections.axon_projections[axis];

            if ((position[axis] >= axis_limits[0] - distance_to_be_inside) && 
                (position[axis] <= axis_limits[1] + distance_to_be_inside))
            {
                ++count_isnear;
            }
        }
        if (count_isnear == 3)
        {
            return tuple<string, int>{"dendrite", i};
        }
    }

    return tuple<string, int>{"none", -1};
}


bool Neuron::checkCollision(Walker &walker, Eigen::Vector3d &step_dir, double &step_lenght, Collision &colision)
{
    // Check if collision with soma
    if(soma.checkCollision(walker, step_dir, step_lenght, colision))
    {
        return true;
    } 
    // Check if collision with dendrites
    for(int i=0 ; i < nb_dendrites ; ++i)
    {
        if(dendrites[i].checkCollision(walker, step_dir, step_lenght, colision))
        {
            return true;
        }  
    } 
    return false;
}


bool Neuron::intersection_sphere_vector(double &intercept1, double &intercept2, Dynamic_Sphere &sphere, Vector3d &step_dir, double &step_length, Vector3d &traj_origin, double &c){
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // A sphere and a line can intersect in 3 ways : no intersect, 1 intersect (tangent), 2 intersect

    Vector3d m = traj_origin - sphere.center;
    double rad = sphere.radius;

    // collision distance
    double d_ = m.norm() - rad;

    // If the minimum distance from the walker to the sphere is more than
    // the actual step size, we can discard this collision.
    if(d_> EPS_VAL){
        if(d_ > step_length+barrier_tickness){
            return false;
        }
    }

    double a = 1;
    double b = (m.dot(step_dir));
    c = m.dot(m) - rad*rad;

    
    double discr = b*b - a*c;

    if (discr < 0.0 ){
        return false;
    }

    intercept1 = (-b + sqrt(discr))/(a);
    intercept2 = (-b - sqrt(discr))/(a);

    return true;

}

tuple<string, int, int> Neuron::closest_sphere_dichotomy(Walker &walker, double &step_lenght, double barrier_thickness)
{
    Vector3d O;
    walker.getVoxelPosition(O);

    string neuron_part; // "soma", "dendrite" or "none"
    int part_id; // id of the soma or dendrite. -1 if not in neuron
    tie(neuron_part, part_id) = isNearNeuron(O, barrier_thickness);

    if (neuron_part=="dendrite")
    {
        int number_spheres = dendrites[part_id].spheres.size();
        int i=0;
        int i_last = number_spheres-1;
        bool stop = false;
        int count= 0;
        int half_way;
        double first_distance;
        double last_distance;
        while ((i_last-i)>1)
        {
            half_way       = int((i_last+i)/2);
            first_distance = (dendrites[part_id].spheres[i].center-O).norm() - dendrites[part_id].spheres[i].radius;
            last_distance  = (dendrites[part_id].spheres[i_last].center-O).norm() - dendrites[part_id].spheres[i_last].radius;
            if(first_distance < last_distance)
            {
                i_last = half_way;
            } 
            else
            {
                i = half_way;
            } 
            count += 1;
        }
        if (first_distance < last_distance)
        {
            return {neuron_part, part_id, i};
        }
        else
        {
            return {neuron_part, part_id, i_last};
        }  
    }
    else if (neuron_part == "soma")
    {
        return {neuron_part, part_id, 0};
    }
} 


void Neuron::set_axons(std::vector<Axon> axons_to_add, int axon_id)
{
    if (axons_to_add.size() != 0)
    {
        dendrites = axons_to_add;
    }
}

void Neuron::generateSpanRadius(int lower_bound, int upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist_span_radius(lower_bound, upper_bound);

    // Generate int number in [lb, ub], in [mm]
    span_radius  = dist_span_radius(rng);
    span_radius /= 10;
}