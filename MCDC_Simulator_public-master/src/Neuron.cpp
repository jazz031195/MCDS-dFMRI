#include "Neuron.h"
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

    constexpr uint8_t ub = 25; // upper bound for nb_dendrites
    constexpr uint8_t lb = 15;  // lower bound for nb_dendrites

    uniform_int_distribution<mt19937::result_type> dist_dendrites(lb, ub);
    // Generate int number in [lb, ub]
    nb_dendrites = 1;//dist_dendrites(rng);

    // // Create a random span radius and set its value to this
    // generateSpanRadius();

    dendrites.clear();

}

Neuron::~Neuron()
{
    nb_neurons--;
}

Neuron::Neuron(vector<Dendrite> const& dendrites_, Dynamic_Sphere const& soma_) : Neuron()
{
    dendrites    = dendrites_;     
    assert(dendrites_.size() < std::numeric_limits<decltype(nb_dendrites)>::max());
    nb_dendrites = dendrites_.size();
    soma         = soma_;
}

Neuron::Neuron(Vector3d const& soma_center, double const& soma_radius=5e-3) : Neuron()
{
    soma = Dynamic_Sphere(soma_center, soma_radius);
}

Neuron::Neuron(vector<Dendrite> const& dendrites_, Vector3d const& soma_center, double const& soma_radius=5e-3) : Neuron()
{
    dendrites    = dendrites_; 
    nb_dendrites = dendrites_.size();        
    soma         = Dynamic_Sphere(soma_center, soma_radius);
}

Neuron::Neuron(Neuron const& neuron) : nb_dendrites(neuron.nb_dendrites), span_radius(neuron.span_radius), dendrites(neuron.dendrites), soma(neuron.soma) 
{
    id = nb_neurons++;
}

double Neuron::minDistance(Walker const& walker) const
{
    vector<double> distances  = Distances_to_Spheres(walker);
    double min = *min_element(begin(distances), end(distances));

    return min;
}

vector<double> Neuron::Distances_to_Spheres(Walker const& w) const
{

    Vector3d O;
    w.getVoxelPosition(O);
    return Distances_to_Spheres(O);
}

double Neuron::minDistance(Eigen::Vector3d const& pos) const
{

    vector<double> distances = Distances_to_Spheres(pos);
    double min = *min_element(begin(distances), end(distances));
    return min;
}

vector<double> Neuron::Distances_to_Spheres(Vector3d const& pos) const
{
 
    // First check distance to soma
    Vector3d m = pos - soma.center;
    double distance_to_sphere = m.norm() - soma.radius;
    vector<double> distances {distance_to_sphere};

    // Then iterate through dendrites
    for (uint8_t i=0; i < dendrites.size(); ++i)
    {
        // Then iterate through the subbranches of each dendrite
        for (size_t j=0; j < dendrites[i].subbranches.size(); j++)
        {   
            // Then iterate through the spheres of a subbranch
            for (size_t k=0; k < dendrites[i].subbranches[j].spheres.size(); k++)
            {
                Vector3d m = pos - dendrites[i].subbranches[j].spheres[k].center;
                double distance_to_sphere = m.norm() - dendrites[i].subbranches[j].spheres[k].radius;
                distances.push_back(distance_to_sphere);
            }    
        }
    }
    return distances;
}

bool Neuron::isPosInsideNeuron(Eigen::Vector3d const& position, double const& distance_to_be_inside, bool const& swell_, int& in_soma_index, int& in_dendrite_index, int& in_subbranch_index)
{
    // if position is in box with Dendrite inside
    string neuron_part; // "dendrite" or "none"
    vector<int> part_id;        // id of the dendrite. -1 if not in neuron
    // distance_to_be_inside = position.sphere.radius + constant
    tie(neuron_part, part_id) = isNearDendrite(position, distance_to_be_inside);
    // cout << neuron_part << endl;
    std::vector<int> sphere_ids;
    if (!(neuron_part == "none"))
    {
        for (size_t p=0; p < part_id.size(); p++)
            for (size_t b=0; b < dendrites[p].subbranches.size(); b++)
            {
                if (dendrites[p].subbranches[b].isPosInsideAxon(position, distance_to_be_inside, false, sphere_ids))
                {
                    in_dendrite_index  = p;
                    in_subbranch_index = b;
                    in_soma_index      = -1;
                    // cout << "in dendrite" << endl;
                    return true;
                }
            }
        
    }

    int soma_id;
    // "soma" or "none". 0 if in soma, -1 otherwise
    tie(neuron_part, soma_id) = isNearSoma(position, distance_to_be_inside);
    if (neuron_part == "soma")
    {
        if (soma.isInside(position, distance_to_be_inside))
        {
            in_soma_index      = 0;
            in_dendrite_index  = -1;
            in_subbranch_index = -1;
            // cout << "in soma" << endl;
            return true;
        }
    }
    in_dendrite_index  = -1;
    in_subbranch_index = -1;
    in_soma_index      = -1;
    return false;
} 

tuple<string, int> Neuron::isNearSoma(Vector3d const& position, double const& distance_to_be_inside) const
{
    int count_isnear = 0;
    // Check soma box
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
        if ((position[axis] >= soma.center[axis] - soma.radius ) &&
            (position[axis] <= soma.center[axis] + soma.radius ))
        {
            ++count_isnear;
        }
    }
    if (count_isnear == 3)
        return tuple<string, int>{"soma", soma.id};
    else
        return tuple<string, int>{"none", -1};
}

tuple<string, vector<int>> Neuron::isNearDendrite(Vector3d const& position, double const& distance_to_be_inside) const
{
    // TODO [ines] : it may be possible that it is near 2 dendrites, and that the first one checked will be returned, even though
    // the walker is in the other...

    // Check each dendrite's box
    int count_isnear = 0;
    vector <int> dendrite_ids;
    for (unsigned int i = 0; i < dendrites.size(); ++i)
    {
        count_isnear = 0;
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
            Vector2d axis_limits = dendrites[i].projections.axon_projections[axis];
            if ((position[axis] >= axis_limits[0] ) &&
                (position[axis] <= axis_limits[1] ))
            {
                ++count_isnear;
            }

        }
        // Inside the box around dendrite
        if (count_isnear == 3)
            dendrite_ids.push_back(i);
            

    }
    if (dendrite_ids.size() > 0)
        return tuple<string, vector<int>>{"dendrite", dendrite_ids};

    return tuple<string, vector<int>>{"none", {-1}};
}

bool Neuron::checkCollision(Walker &walker, Vector3d const& step_dir, double const& step_lenght, Collision &colision)
{
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d next_step = step_dir * step_lenght + O;
    vector<int> sphere_ids;
    
    // If inside soma
    if (walker.in_soma_index == 0)
    {        
        // But next step isn't
        if (!soma.isInside(next_step, -barrier_tickness))
        {
            // Find in which dendrite it can go
            int closest_dendrite_id = closest_dendrite_from_soma(next_step, step_lenght);

            // If inside closest_dendrite_id dendrite => no collision
            // TODO : check if enough to check only the first subbranch [ines]
            if ((closest_dendrite_id >= 0) &&
                dendrites[closest_dendrite_id].subbranches[0].isPosInsideAxon(next_step, -barrier_tickness, false, sphere_ids))
            {
                colision.type = Collision::null;
                walker.in_dendrite_index = closest_dendrite_id;
                walker.in_subbranch_index= 0;
                walker.in_soma_index     = -1;
                cout << "next step in dendrite" << endl;
                return false;
            }
            // Or in extra => collision
            else if (soma.checkCollision(walker, step_dir, step_lenght, colision))
            {
                // cout << "exit soma " << endl;
                return true;
            }
                
        }
        // If the next step is also in the soma => no collision
        else
        {
            // cout << "stays in soma" << endl;
            walker.in_soma_index      = 0;
            walker.in_dendrite_index  = -1;
            walker.in_subbranch_index = -1;
            colision.type = Collision::null;
            return false;
        }
    }
    // Check if collision with soma. Can be internal or external collisions.
    if (walker.location == Walker::extra)
        if (soma.checkCollision(walker, step_dir, step_lenght, colision))
        {
            return true;
        }

    // If inside a dendrite
    if (walker.in_dendrite_index >= 0)
    {   
        // cout << "in dendrite " << endl;
        vector<int> branching_id = closest_subbranch(next_step, walker.in_dendrite_index, walker.in_subbranch_index, step_lenght);
        const auto& subbranches  = dendrites[walker.in_dendrite_index].subbranches;
        bool next_step_in_subbranch = subbranches[walker.in_subbranch_index].isPosInsideAxon(next_step, -barrier_tickness, false, sphere_ids);

        // The next step is in the same subbranch => no collision
        if(branching_id[0] == -1 && next_step_in_subbranch)
        {
            // cout << "still in dendrite" << endl;
            colision.type = Collision::null;
            walker.in_soma_index = -1;
            return false;
        }
        // The next step is in the outside => collision
        else if(branching_id[0] == -1 && !next_step_in_subbranch)
        {
            cout << "exit this dendrite" << endl;
            if(subbranches[walker.in_subbranch_index].checkCollision(walker, step_dir, step_lenght, colision))
            {   
                return true;
            }
            cout << "why here" << endl;
        }
        // The walker is close to soma
        else if (branching_id[0] == 0 && branching_id.size() == 1)
        {
            // Next step in the same subbranch
            if (next_step_in_subbranch)
            {
                // cout << "stays in the branch" << endl;
                colision.type = Collision::null;
                walker.in_soma_index = -1;
                return false;
            }
            // Next step is in soma
            else if (soma.isInside(next_step, -barrier_tickness))
            {
                colision.type            = Collision::null;
                walker.in_soma_index     = 0;
                walker.in_dendrite_index = -1;
                walker.in_subbranch_index= -1;
                cout << "next step in soma" << endl;
                return false;
            }
            // If in extra => collision
            else if (subbranches[walker.in_subbranch_index].checkCollision(walker, step_dir, step_lenght, colision, soma))
            {
                cout << "exits the branch" << endl;
                return true;
            }
                
        }
        // Next step is in another subbranch
        else if((branching_id.size() > 1) && (subbranches[branching_id[0]].isPosInsideAxon(next_step, -barrier_tickness, false, sphere_ids)))
        {
            cout << "goes to another branch" << endl;
            colision.type = Collision::null;
            walker.in_soma_index      = -1;
            walker.in_subbranch_index = branching_id[0];
            return false;
        }
        // Next step is in another subbranch
        else if((branching_id.size() > 1) && (subbranches[branching_id[1]].isPosInsideAxon(next_step, -barrier_tickness, false, sphere_ids)))
        {
            cout << "goes to another branch" << endl;
            colision.type = Collision::null;
            walker.in_soma_index      = -1;
            walker.in_subbranch_index = branching_id[1];
            return false;
        }
    }
    // Check if collision with dendrites from outside.
    if (walker.location == Walker::extra)
    {
        for (uint8_t i = 0; i < dendrites.size(); ++i)
        {
            if (dendrites[i].checkCollision(walker, step_dir, step_lenght, colision))
                return true;
        }
    }
    return false;
}

vector <int> Neuron::closest_subbranch(Vector3d const& position, int const& dendrite_id, int const& subbranch_id, double const& step_length)
{
    const auto& subbranch = dendrites[dendrite_id].subbranches[subbranch_id];
    int size_subbranch    = subbranch.spheres.size() - 1;
    double eps = 0.0005;
    double distance_to_proximal_end = (subbranch.spheres[0].center - position).norm() - eps;
    double distance_to_distal_end   = (subbranch.spheres[size_subbranch].center - position).norm() - eps;

    if(distance_to_proximal_end < step_length)
    {
        // The subbranch id starts at 1 but the indices in c++ start at 0
        vector<int> proximal_branching = subbranch.proximal_branching;
        if(proximal_branching.size() > 1)
        {
           proximal_branching[0] -= 1;  
           proximal_branching[1] -= 1;   
           cout << subbranch.proximal_branching[0] << subbranch.proximal_branching[1] << endl;
        }
        return proximal_branching;
    }    
    else if (distance_to_distal_end < step_length)
    {
        // The subbranch id starts at 1 but the indices in c++ start at 0
        vector<int> distal_branching = subbranch.distal_branching;
        
        // No need to check the size, because there are always 2 distal branching
        distal_branching[0] -= 1;  
        distal_branching[1] -= 1; 
        cout << subbranch.distal_branching[0] << subbranch.distal_branching[1] << endl;  
        
        return distal_branching;
    }   
    else
        return {-1};
}

int Neuron::closest_dendrite_from_soma(Vector3d const& position, double const& step_length)
{
    int smaller_dist = 1000;
    int closer_dendrite = -1;
    for (size_t i = 0; i < dendrites.size(); ++i)
    {
        // TODO : there was a segfault here once ... [ines]
        int distance_to_axon = (position - dendrites[i].subbranches[0].spheres[0].center).norm();
        if (distance_to_axon < smaller_dist && distance_to_axon <= step_length)
        {
            smaller_dist    = distance_to_axon;
            closer_dendrite = i;
        }
    }
    return closer_dendrite;
}

bool Neuron::intersection_sphere_vector(double &intercept1, double &intercept2, Dynamic_Sphere const& sphere, Vector3d const& step_dir, double const& step_length, Vector3d const& traj_origin, double &c) const
{
    Vector3d m = traj_origin - sphere.center;
    double rad = sphere.radius;

    // collision distance
    double d_ = m.norm() - rad;

    // If the minimum distance from the walker to the sphere is more than
    // the actual step size, we can discard this collision.
    if (d_ > EPS_VAL)
        if (d_ > step_length + barrier_tickness)
            return false;

    double a = 1;
    double b = m.dot(step_dir);
    c = m.dot(m) - rad * rad;

    double discr = b * b - a * c;

    if (discr < 0.0)
        return false;

    intercept1 = (-b + sqrt(discr)) / (a);
    intercept2 = (-b - sqrt(discr)) / (a);

    return true;
}

// tuple<string, int, int> Neuron::closest_sphere(Walker const& walker, double const& barrier_thickness) const
// {
//     Vector3d walker_pos;
//     walker.getVoxelPosition(walker_pos);

//     string neuron_part; // "soma", "dendrite" or "none"
//     int part_id;        // id of the soma or dendrite. -1 if not in neuron
//     tie(neuron_part, part_id) = isNearNeuron(walker_pos, barrier_thickness, false, false);

//     if (neuron_part == "dendrite")
//     {
//         int i = closest_sphere_dichotomy(part_id, walker_pos);
//         return {neuron_part, part_id, i};
//     }
//     else if (neuron_part == "soma")
//         return {neuron_part, part_id, 0};
//     else
//         return {"none", -1, -1}; // TODO: to check
// }

// int Neuron::closest_sphere_dichotomy(int const& dendrite_id, Vector3d const& walker_pos) const
// {
//     int i = 0;
//     int number_spheres = dendrites[dendrite_id].spheres.size();
//     int i_last = number_spheres - 1;
//     int half_way;
//     double first_distance;
//     double last_distance;
//     while ((i_last - i) > 1)
//     {
//         half_way = int((i_last + i) / 2);
//         first_distance = (dendrites[dendrite_id].spheres[i].center - walker_pos).norm() - dendrites[dendrite_id].spheres[i].radius;
//         last_distance = (dendrites[dendrite_id].spheres[i_last].center - walker_pos).norm() - dendrites[dendrite_id].spheres[i_last].radius;
//         if (first_distance < last_distance)
//             i_last = half_way;
//         else
//             i = half_way;
//     }
//     if (first_distance < last_distance)
//         return i;

//     return i_last;
// }

void Neuron::set_dendrites(std::vector<Dendrite> const& dendrites_to_add)
{
    if (dendrites_to_add.size() != 0)
        dendrites = dendrites_to_add;
}

void Neuron::add_dendrite(Dendrite& dendrite_to_add)
{
    dendrite_to_add.add_projection(dendrite_to_add.id);
    dendrites.push_back(dendrite_to_add);
}

void Neuron::generateSpanRadius(double const& lower_bound, double const& upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<double> dist_span_radius(lower_bound, upper_bound);
    cout << dist_span_radius(rng) << endl;
    // Generate int number in [lb, ub], in [mm]
    span_radius = dist_span_radius(rng);
}

vector<double> Neuron::get_Volume() const
{
    double VolumeSoma = 0;
    double VolumeDendrites = 0;
    // Calculate the volume of the soma
    VolumeSoma += 4/3*M_PI*pow(soma.radius, 3);

    // Calculate the cylindrical volume of each dendrite
    for (uint8_t j = 0; j < nb_dendrites; j++)
    {
        VolumeDendrites += dendrites[j].volumeDendrite();
    } 

    return {VolumeSoma, VolumeDendrites};
}