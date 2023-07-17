#include "doctest.h"

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
    constexpr uint8_t lb = 15; // lower bound for nb_dendrites

    uniform_int_distribution<mt19937::result_type> dist_dendrites(lb, ub);
    // Generate int number in [lb, ub]
    nb_dendrites = 1; // dist_dendrites(rng);

    // // Create a random span radius and set its value to this
    // generateSpanRadius();

    dendrites.clear();
}

Neuron::~Neuron()
{
    nb_neurons--;
}

Neuron::Neuron(vector<Dendrite> const &dendrites_, Dynamic_Sphere const &soma_) : Neuron()
{
    dendrites = dendrites_;
    assert(dendrites_.size() < std::numeric_limits<decltype(nb_dendrites)>::max());
    nb_dendrites = dendrites_.size();
    soma = soma_;
}

Neuron::Neuron(Vector3d const &soma_center, double const &soma_radius, int const &neuron_id) : Neuron()
{
    soma = Dynamic_Sphere(soma_center, soma_radius, neuron_id);
}

Neuron::Neuron(vector<Dendrite> const &dendrites_, Vector3d const &soma_center, double const &soma_radius = 10e-3) : Neuron()
{
    dendrites = dendrites_;
    nb_dendrites = dendrites_.size();
    soma = Dynamic_Sphere(soma_center, soma_radius, id);
}

Neuron::Neuron(Neuron const &neuron) : nb_dendrites(neuron.nb_dendrites), span_radius(neuron.span_radius), dendrites(neuron.dendrites), soma(neuron.soma)
{
    id = nb_neurons++;
}
double Neuron::minDistance(Walker const &walker) const
{
    vector<double> distances = Distances_to_Spheres(walker);
    double min = *min_element(begin(distances), end(distances));

    return min;
}

vector<double> Neuron::Distances_to_Spheres(Walker const &w) const
{

    Vector3d O;
    w.getVoxelPosition(O);
    return Distances_to_Spheres(O);
}

double Neuron::minDistance(Eigen::Vector3d const &pos) const
{

    vector<double> distances = Distances_to_Spheres(pos);
    double min = *min_element(begin(distances), end(distances));
    return min;
}

vector<double> Neuron::Distances_to_Spheres(Vector3d const &pos) const
{

    // First check distance to soma
    Vector3d m = pos - soma.center;
    double distance_to_sphere = m.norm() - soma.radius;
    vector<double> distances{distance_to_sphere};

    // Then iterate through dendrites
    for (uint8_t i = 0; i < dendrites.size(); ++i)
    {
        // Then iterate through the subbranches of each dendrite
        for (size_t j = 0; j < dendrites[i].subbranches.size(); j++)
        {
            // Then iterate through the spheres of a subbranch
            for (size_t k = 0; k < dendrites[i].subbranches[j].spheres.size(); k++)
            {
                Vector3d m = pos - dendrites[i].subbranches[j].spheres[k].center;
                double distance_to_sphere = abs(m.norm() - dendrites[i].subbranches[j].spheres[k].radius);
                distances.push_back(distance_to_sphere);
            }
        }
    }
    return distances;
}

TEST_CASE("check_isPosInsideNeuron")
{
    Vector3d center(0.05, 0.05, 0.05);
    double radius_soma = 10e-3;
    double radius_dendrite = 0.5e-3;
    int sphere_id = 0;

    Neuron neuron(center, radius_soma, 0);
    Dendrite dendrite;

    int branch_id = 1;
    vector<Dynamic_Sphere> spheres_list;
    spheres_list.push_back(Dynamic_Sphere(center, radius_soma, sphere_id));
    for (size_t i = 0; i < 10; ++i)
    {
        Vector3d next_center(center[0] + radius_soma + i * radius_dendrite / 4, center[1], center[2]);
        Dynamic_Sphere sphere_to_add(next_center, radius_dendrite, 0, false, branch_id, i, 1);
        spheres_list.push_back(sphere_to_add);
    }

    Vector3d last_center(spheres_list[spheres_list.size()-1].center);
    Vector3d begin;
    vector <int> proximal_end = {1};
    vector <int> distal_end   = {2, 3};
    Axon subbranch(branch_id, radius_dendrite, begin, begin, 0, false, 1, proximal_end, distal_end);
    subbranch.set_spheres(spheres_list);
    dendrite.add_subbranch(subbranch);

    // Add branching
    Vector3d branching_direction(1, 1, 0);
    branching_direction = branching_direction.normalized();

    for(size_t rep=0; rep < 2; rep++)
    {
        spheres_list.clear();
        for (size_t i = 0; i < 10; ++i)
        {
            Vector3d next_center(last_center[0] + i * radius_dendrite / 4 * branching_direction[0],
                                 last_center[1] + i * radius_dendrite / 4 * branching_direction[1],
                                 last_center[2] + i * radius_dendrite / 4 * branching_direction[2]);
            Dynamic_Sphere sphere_to_add(next_center, radius_dendrite, 0, false, branch_id + rep + 1, i, 1);
            spheres_list.push_back(sphere_to_add);
        }
        vector<int> proximal_end {1, int(3-rep)};
        vector<int> distal_end {int(4+rep), int(5+rep)};
        Axon subbranch_(branch_id + rep + 1, radius_dendrite, begin, begin, 0, false, 1, proximal_end, distal_end);
        subbranch_.set_spheres(spheres_list);
        dendrite.add_subbranch(subbranch_);
        branching_direction = - branching_direction;
    }
    neuron.add_dendrite(dendrite);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);
    int dummy;
    vector<int> dummy1;
    SUBCASE("soma")
    {
        Vector3d somaCenter = neuron.soma.center;
        double somaRadius   = neuron.soma.radius;
        
        for(size_t nb_t=0; nb_t < 10; nb_t++)
        {
            // Check value inside soma
            double probaRadius = double(udist(gen));
            
            SUBCASE("inside soma")
            {      
                double theta = 2 * M_PI * udist(gen);
                double phi = acos(1 - 2 * udist(gen));
                double x = sin(phi) * cos(theta) * probaRadius * somaRadius + somaCenter[0];
                double y = sin(phi) * sin(theta) * probaRadius * somaRadius + somaCenter[1];
                double z = cos(phi) * probaRadius * somaRadius + somaCenter[2];
                Vector3d pos_temp = {x,y,z};

                bool inSoma_formula = pow(x-somaCenter[0], 2) + pow(y-somaCenter[1], 2) + pow(z-somaCenter[2], 2) <= somaRadius;
                bool inSoma_function= neuron.isPosInsideNeuron(pos_temp, EPS_VAL, false, dummy, dummy, dummy, dummy1);
                CHECK_EQ(inSoma_formula, inSoma_function);
            }
            SUBCASE("at soma surface")
            {
                double theta = 2 * M_PI * udist(gen);
                double phi = acos(1 - 2 * udist(gen));
                double x = sin(phi) * cos(theta) * somaRadius + somaCenter[0];
                double y = sin(phi) * sin(theta) * somaRadius + somaCenter[1];
                double z = cos(phi) * somaRadius + somaCenter[2];
                Vector3d pos_temp = {x,y,z};

                bool inSoma_formula = pow(x-somaCenter[0], 2) + pow(y-somaCenter[1], 2) + pow(z-somaCenter[2], 2) <= somaRadius;
                bool inSoma_function= neuron.isPosInsideNeuron(pos_temp, EPS_VAL, false, dummy, dummy, dummy, dummy1);
                CHECK_EQ(inSoma_formula, inSoma_function);
            }
        }     
    }
    SUBCASE("in dendrite")
    {
        for(size_t d=0; d < neuron.dendrites.size(); ++ d)
        {
            for(size_t b=0; b < neuron.dendrites[d].subbranches.size(); ++ b)
            {
                auto subbranch = neuron.dendrites[d].subbranches[b];
                for(size_t s=0; s < subbranch.spheres.size(); ++ s)
                {
                    Vector3d sphereCenter = subbranch.spheres[s].center;
                    double sphereRadius   = subbranch.spheres[s].radius;
                    SUBCASE("inside dendrite")
                    {
                        // Check value inside spheres
                        double probaRadius = double(udist(gen));
                        
                        double theta = 2 * M_PI * udist(gen);
                        double phi = acos(1 - 2 * udist(gen));
                        double x = sin(phi) * cos(theta) * probaRadius * sphereRadius + sphereCenter[0];
                        double y = sin(phi) * sin(theta) * probaRadius * sphereRadius + sphereCenter[1];
                        double z = cos(phi) * probaRadius * sphereRadius + sphereCenter[2];
                        Vector3d pos_temp = {x,y,z};

                        bool inSphere_formula = pow(x-sphereCenter[0], 2) + pow(y-sphereCenter[1], 2) + pow(z-sphereCenter[2], 2) <= sphereRadius;
                        bool inSphere_function= neuron.isPosInsideNeuron(pos_temp, EPS_VAL, false, dummy, dummy, dummy, dummy1);
                        CHECK_EQ(inSphere_formula, inSphere_function);
                    }
                    SUBCASE("at dendrite surface")
                    {
                        double theta = 2 * M_PI * udist(gen);
                        double phi = acos(1 - 2 * udist(gen));
                        double x = sin(phi) * cos(theta)  * sphereRadius + sphereCenter[0];
                        double y = sin(phi) * sin(theta)  * sphereRadius + sphereCenter[1];
                        double z = cos(phi) * sphereRadius + sphereCenter[2];
                        Vector3d pos_temp = {x,y,z};

                        bool inSphere_formula = pow(x-sphereCenter[0], 2) + pow(y-sphereCenter[1], 2) + pow(z-sphereCenter[2], 2) <= sphereRadius;
                        bool inSphere_function= neuron.isPosInsideNeuron(pos_temp, EPS_VAL, false, dummy, dummy, dummy, dummy1);
                        CHECK_EQ(inSphere_formula, inSphere_function);
                    }
                }
            }
        }
    }
}

bool Neuron::isPosInsideNeuron(Eigen::Vector3d const &position, double const &distance_to_be_inside, bool const &swell_, int &in_soma_index, int &in_dendrite_index, int &in_subbranch_index, vector<int> &in_sph_index)
{
    vector<int> sphere_ids;
    for (size_t p = 0; p < dendrites.size(); p++)
    {
        for (size_t b = 0; b < dendrites[p].subbranches.size(); b++)
        {
            for (size_t s = 0; s < dendrites[p].subbranches[b].spheres.size(); s++)
            {
                Vector3d somaCenter = dendrites[p].subbranches[b].spheres[s].center;
                double somaRadius = dendrites[p].subbranches[b].spheres[s].radius;

                if (pow(position[0] - somaCenter[0], 2) + pow(position[1] - somaCenter[1], 2) + pow(position[2] - somaCenter[2], 2) <= somaRadius + distance_to_be_inside)
                {
                    in_dendrite_index = p;
                    in_subbranch_index = b;
                    in_soma_index = -1;
                    in_sph_index = sphere_ids;
                    return true;
                }
            }
        }
    }

    if (soma.isInside(position, distance_to_be_inside))
    {
        in_soma_index = 0;
        in_dendrite_index = -1;
        in_subbranch_index = -1;
        in_sph_index.clear();
        return true;
    }

    // // if position is in box with Dendrite inside
    // string neuron_part; // "dendrite" or "none"
    // vector<int> part_id;        // id of the dendrite. -1 if not in neuron
    // // distance_to_be_inside = position.sphere.radius + constant
    // tie(neuron_part, part_id) = isNearDendrite(position, distance_to_be_inside);
    // // cout << neuron_part << endl;
    // std::vector<int> sphere_ids;
    // if (!(neuron_part == "none"))
    // {
    //     for (size_t p=0; p < part_id.size(); p++)
    //     {
    //         for (size_t b=0; b < dendrites[p].subbranches.size(); b++)
    //         {
    //             if (dendrites[p].subbranches[b].isPosInsideAxon(position, distance_to_be_inside, sphere_ids))
    //             {
    //                 in_dendrite_index  = p;
    //                 in_subbranch_index = b;
    //                 in_soma_index      = -1;
    //                 in_sph_index       = sphere_ids;
    //                 // cout << "in dendrite" << endl;
    //                 return true;
    //             }
    //         }
    //     }
    // }

    // int soma_id;
    // // "soma" or "none". 0 if in soma, -1 otherwise
    // tie(neuron_part, soma_id) = isNearSoma(position, distance_to_be_inside);
    // if (neuron_part == "soma")
    // {
    //     if (soma.isInside(position, distance_to_be_inside))
    //     {
    //         in_soma_index      = 0;
    //         in_dendrite_index  = -1;
    //         in_subbranch_index = -1;
    //         in_sph_index.clear();
    //         // cout << "in soma" << endl;
    //         return true;
    //     }
    // }
    // in_dendrite_index  = -1;
    // in_subbranch_index = -1;
    // in_soma_index      = -1;
    // in_sph_index.clear();
    return false;
}

tuple<string, int> Neuron::isNearSoma(Vector3d const &position, double const &distance_to_be_inside) const
{
    int count_isnear = 0;
    // Check soma box
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
        if ((position[axis] >= soma.center[axis] - soma.radius) &&
            (position[axis] <= soma.center[axis] + soma.radius))
        {
            ++count_isnear;
        }
    }
    if (count_isnear == 3)
        return tuple<string, int>{"soma", soma.id};
    else
        return tuple<string, int>{"none", -1};
}

tuple<string, vector<int>> Neuron::isNearDendrite(Vector3d const &position, double const &distance_to_be_inside) const
{
    // TODO [ines] : it may be possible that it is near 2 dendrites, and that the first one checked will be returned, even though
    // the walker is in the other...

    // Check each dendrite's box
    int count_isnear = 0;
    vector<int> dendrite_ids;
    for (unsigned int i = 0; i < dendrites.size(); ++i)
    {
        count_isnear = 0;
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
            Vector2d axis_limits = dendrites[i].projections.axon_projections[axis];
            if ((position[axis] >= axis_limits[0]) &&
                (position[axis] <= axis_limits[1]))
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

vector<Dynamic_Sphere *> Neuron::find_neighbor_spheres(Walker &walker, Vector3d const &next_step, double const &step_length)
{
    vector<Dynamic_Sphere *> spheres_list;
    spheres_list.clear();

    spheres_list.push_back(new Dynamic_Sphere(soma));

    if (walker.in_soma_index == 0)
    {
        cout << "in_soma_index == 0" << endl;
        return spheres_list;
    }
    if ((walker.in_dendrite_index == -1) || (walker.in_subbranch_index == -1))
    {
        cout << "weird" << endl;
        return spheres_list;
    }

    const auto subbranch = dendrites[walker.in_dendrite_index].subbranches[walker.in_subbranch_index];
    for (size_t i = 0; i < subbranch.spheres.size(); i++)
    {
        spheres_list.push_back(new Dynamic_Sphere(subbranch.spheres[i]));
    }

    // // We are in the soma
    // if(walker.in_soma_index == 0)
    // {
    //     spheres_list.push_back(&soma);
    //    // Find in which dendrite it can go
    //     int closest_dendrite_id = closest_dendrite_from_soma(next_step, step_length);

    //     // If inside closest_dendrite_id dendrite => no collision
    //     // TODO : check if enough to check only the first subbranch [ines]
    //     if (closest_dendrite_id >= 0)
    //     {
    //         for(size_t i=0; i < 5; ++i)
    //         {
    //             spheres_list.push_back(&dendrites[closest_dendrite_id].subbranches[0].spheres[i]);
    //         }
    //     }
    // }
    // if(walker.in_dendrite_index >= 0)
    // {
    //     vector<int> branching_id = closest_subbranch(next_step, walker.in_dendrite_index, walker.in_subbranch_index, step_length);

    //     // The walker is close to soma
    //     if (branching_id[0] == 0 && branching_id.size() == 1)
    //     {
    //         cout << "close to soma " << endl;
    //         spheres_list.push_back(&soma);
    //         for(size_t i=0; i < 5; ++i)
    //         {
    //             spheres_list.push_back(&dendrites[walker.in_dendrite_index].subbranches[0].spheres[i]);
    //         }

    //     }
    //     // We are in a middle of a subbranch
    //     else if (branching_id[0] == -1)
    //     {
    //         cout << "middle" << endl;
    //         int middle_idx = walker.in_sph_index.size()/2;
    //         // Add the current sphere
    //         spheres_list.push_back(&dendrites[walker.in_dendrite_index].subbranches[0].spheres[middle_idx]);

    //         for(size_t i=0; i < 5; ++i)
    //         {
    //             // cout << middle_idx - i << endl;
    //             // cout << middle_idx + i << endl;
    //             // cout << dendrites[walker.in_dendrite_index].subbranches[0].spheres.size() << endl;
    //             // Add spheres before the current
    //             if(middle_idx >= i)
    //             {
    //                 spheres_list.push_back(&dendrites[walker.in_dendrite_index].subbranches[0].spheres[middle_idx - i]);
    //             }
    //             // Add spheres after the current
    //             if((middle_idx + i) < dendrites[walker.in_dendrite_index].subbranches[0].spheres.size())
    //                 spheres_list.push_back(&dendrites[walker.in_dendrite_index].subbranches[0].spheres[middle_idx + i]);
    //         }
    //     }
    //     // Distal end of the neuron reached
    //     else if (branching_id[0] == -2)
    //     {

    //         cout << "distal end" << walker.in_subbranch_index << endl;
    //         auto subbranch = dendrites[walker.in_dendrite_index].subbranches[walker.in_subbranch_index];
    //         for(size_t i=0; i < 5; ++i)
    //         {
    //             spheres_list.push_back(&subbranch.spheres[subbranch.spheres.size()-1-i]);
    //         }
    //     }
    //     // Next step is in another subbranch
    //     else if(branching_id.size() > 1)
    //     {
    //         cout << "another subbranch" << endl;
    //         // distal branching
    //         if(branching_id[0] > walker.in_subbranch_index)
    //         {
    //             cout << "distal" << endl;
    //             cout << walker.in_subbranch_index << endl;
    //             auto current_branch = dendrites[walker.in_dendrite_index].subbranches[walker.in_subbranch_index];
    //             auto child1_branch  = dendrites[walker.in_dendrite_index].subbranches[branching_id[0]];
    //             auto child2_branch  = dendrites[walker.in_dendrite_index].subbranches[branching_id[1]];
    //             for(size_t i=0; i < 5; ++i)
    //             {
    //                 spheres_list.push_back(&current_branch.spheres[current_branch.spheres.size()-1-i]);
    //                 spheres_list.push_back(&child1_branch.spheres[i]);
    //                 spheres_list.push_back(&child2_branch.spheres[i]);
    //             }
    //         }
    //         // proximal branching
    //         else
    //         {
    //             cout << "proximal" << endl;
    //             cout << walker.in_subbranch_index << endl;
    //             auto current_branch = dendrites[walker.in_dendrite_index].subbranches[walker.in_subbranch_index];
    //             auto parent_branch  = dendrites[walker.in_dendrite_index].subbranches[branching_id[0]];
    //             auto child_branch   = dendrites[walker.in_dendrite_index].subbranches[branching_id[1]];
    //             for(size_t i=0; i < 5; ++i)
    //             {
    //                 spheres_list.push_back(&current_branch.spheres[i]);
    //                 spheres_list.push_back(&parent_branch.spheres[parent_branch.spheres.size()-1-i]);
    //                 spheres_list.push_back(&child_branch.spheres[i]);
    //             }
    //         }
    //     }

    // }

    return spheres_list;
}

bool Neuron::checkCollision(Walker &walker, Vector3d const &step_dir, double const &step_lenght, Collision &colision)
{
    bool isColliding = false;
    // cout << "check" << endl;
    // Vector3d O;
    // walker.getVoxelPosition(O);
    // Vector3d next_step = step_dir;//step_dir * step_lenght + O;
    // vector<int> sphere_ids;

    vector<Dynamic_Sphere *> spheres_list = find_neighbor_spheres(walker, step_dir, step_lenght);
    // cout << spheres_list.size() << endl;
    // cout << walker.pos_r[0] << " " << walker.pos_r[1] << " " << walker.pos_r[2] << endl;

    if (spheres_list.size() > 1)
    {
        isColliding = checkCollision_branching(walker, spheres_list, step_dir, step_lenght, colision);
        for (size_t i = 0; i < spheres_list.size(); i++)
            delete spheres_list[i];
        spheres_list.clear();
    }
    else
        isColliding = soma.checkCollision(walker, step_dir, step_lenght, colision);

    cout << walker.pos_v[0] << " " << walker.pos_v[1] << " " << walker.pos_v[2] << endl;
    if (!isPosInsideNeuron(walker.pos_v, EPS_VAL, false, walker.in_soma_index, walker.in_dendrite_index, walker.in_subbranch_index, walker.in_sph_index))
    {
        walker.location = Walker::extra;
        cout << "extra" << endl;
    }

    // TODO [ines] : do the extracellular collisions
    // if (walker.location == Walker::extra)
    // {
    //     // Check if collision with soma. Can be internal or external collisions.
    //     if (soma.checkCollision(walker, step_dir, step_lenght, colision))
    //     {
    //         return true;
    //     }
    //     // Check if collision with dendrites from outside.
    //     for (uint8_t i = 0; i < dendrites.size(); ++i)
    //     {
    //         if (dendrites[i].checkCollision(walker, step_dir, step_lenght, colision))
    //             return true;
    //     }
    // }
    if (isColliding)
    {
        cout << "is colliding" << endl;
        return true;
    }
    

    return false;
}

bool check_inside_(vector<double> const &list_c)
{
    for (unsigned i = 0; i < list_c.size(); ++i)
    {
        if (list_c[i] < 1e-10)
            return true;
    }
    return false;
}

bool check_outside_(vector<double> list_c)
{
    for (unsigned i = 0; i < list_c.size(); ++i)
    {
        if (list_c[i] < -1e-10)
        {
            return false;
        }
    }
    return true;
}

TEST_CASE("checkCollision_branching")
{
    Vector3d center(0.05, 0.05, 0.05);
    double radius_soma = 10e-3;
    double radius_dendrite = 0.5e-3;
    int sphere_id = 0;

    vector<Dynamic_Sphere *> spheres_list;
    spheres_list.push_back(new Dynamic_Sphere(center, radius_soma, sphere_id));
    for (size_t i = 0; i < 10; ++i)
    {
        Vector3d next_center(center[0] + radius_soma + i * radius_dendrite / 4, center[1], center[2]);
        spheres_list.push_back(new Dynamic_Sphere(next_center, radius_dendrite, i));
    }
    Collision colision;
    colision.type = Collision::null;
    colision.t = INFINITY_VALUE;

    SUBCASE("no collision")
    {
        Vector3d direction(1, 0, 0);
        double step_length = 5e-4;

        Walker w;
        // w.setNextDirection(direction);
        w.setInitialPosition(center[0] + radius_soma, center[1], center[2]);
        w.location = w.initial_location = Walker::intra;

        bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
        CHECK(!collided);
        CHECK_EQ(colision.type, Collision::null);
        CHECK_EQ(colision.t, INFINITY_VALUE);
    }
    SUBCASE("trivial bounce cycle")
    {
        Vector3d direction(0, 1, 0);
        double step_length = 6e-4;

        Walker w;
        // w.setNextDirection(direction);
        w.setInitialPosition(center[0] + radius_soma, center[1], center[2]);
        w.location = w.initial_location = Walker::intra;

        SUBCASE("trivial hit")
        {
            bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
            CHECK(collided);
            CHECK_EQ(colision.type, Collision::hit);
            CHECK_EQ(colision.t, doctest::Approx(radius_dendrite));
            CHECK_EQ(colision.bounced_direction, -direction);
            CHECK_EQ(colision.colision_point, w.pos_v + radius_dendrite * direction);
            CHECK_EQ(w.status, Walker::free);
        }
        SUBCASE("trivial bounce")
        {
            w.status = Walker::bouncing;
            w.setRealPosition(w.pos_r + radius_dendrite * direction);
            w.setVoxelPosition(w.pos_v + radius_dendrite * direction);

            colision.t = radius_dendrite;
            colision.bounced_direction = -direction;
            colision.colision_point = w.pos_v;

            bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length - radius_dendrite, colision);

            CHECK(!collided);
            CHECK_EQ(colision.type, Collision::null);
            CHECK_EQ(colision.t, doctest::Approx(radius_dendrite));
            CHECK_EQ(colision.bounced_direction, -direction);
            CHECK_EQ(colision.colision_point, w.pos_v);
            CHECK_EQ(w.status, Walker::bouncing);
        }
    }
    // TODO: [ines] more complex cases
    for (size_t i = 0; i < spheres_list.size(); i++)
        delete spheres_list[i];
    spheres_list.clear();
}

bool Neuron::checkCollision_branching(Walker &walker, vector<Dynamic_Sphere *> &spheres_list, Vector3d const &step, double const &step_lenght, Collision &colision)
{

    string message;
    Vector3d O;
    walker.getVoxelPosition(O);

    // distances to intersections
    std::vector<double> dist_intersections;
    // values indicating whether the walker is inside or outside a sphere
    std::vector<double> cs;

    std::vector<int> sph_ids;
    std::vector<double> all_cs;

    for (size_t i = 0; i < spheres_list.size(); ++i)
    {

        // cout << "Sphere : " << i << ", sph_id :" << spheres[i].id <<", Ax :" << id << endl;

        // distances to collision
        double t1;
        double t2;
        double c;
        bool intersect = intersection_sphere_vector(t1, t2, *spheres_list[i], step, step_lenght, O, c);
        // TODO: [ines] check that opposite signs

        double limit_length = -EPS_VAL;

        if (intersect)
        {
            all_cs.push_back(c);
            // // check if the walker has previsouly collided with this sphere
            // // bool same_colliding_object = false;
            // std::vector<int> last_collision;
            // walker.getLastCollision(last_collision);

            // if ((last_collision).size() == 2) {
            // cout << "Last collision of walker : sph : " << walker.last_collision[0] << ", ax : " << walker.last_collision[1] << endl;

            //    if ((i == last_collision[0])&&(id == last_collision[1])) {
            //        same_colliding_object = true;
            // cout << "Same collision as previous !" << endl;
            //    }
            //}
            // cout << "Intersects with sphere " << j <<", t1 :" << t1 << ", t2 : " << t2 <<", c :" << c << endl;
            // if the new position is at edge of axon, at the edge of sphere[i] but inside the neighbours
            bool condition = false;
            if (i == 0)
            {
                condition = !(spheres_list[i + 1]->isInside(step * t1 + O, limit_length));
            }
            else if (i == spheres_list.size() - 1)
            {
                condition = !(spheres_list[i - 1]->isInside(step * t1 + O, limit_length));
            }
            else
            {
                condition = !(spheres_list[i + 1]->isInside(step * t1 + O, limit_length) || spheres_list[i - 1]->isInside(step * t1 + O, limit_length));
                // cout <<  "spheres_list[" << i+1 << "].minDistance(step*t1+O) : "<<spheres_list[i+1].minDistance(step*t1+O) << endl;
                // cout <<  "spheres_list[" << i-1 << "].minDistance(step*t1+O) : "<<spheres_list[i-1].minDistance(step*t1+O) << endl;
                // cout << "this sphere : spheres_list[" << i << "].minDistance(step*t1+O) : "<<spheres_list[i].minDistance(step*t1+O) << endl;
                // cout << "condition : " << condition << endl;
            }
            // condition = true;
            if (condition)
            {

                // if the collision are too close or negative.
                if (Walker::bouncing)
                {
                    if (t1 >= EPS_VAL && t1 <= step_lenght + barrier_tickness)
                    {
                        // cout << "   intersection, sphere (bouncing):" << i << endl;
                        // cout << "       c :" << c << endl;
                        // cout << "       t :" << t1 << endl;
                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
                else
                {
                    if (t1 >= 0 && t1 <= step_lenght + barrier_tickness)
                    {
                        // cout << "   intersection, sphere (bouncing):" << i << endl;
                        // cout << "       c :" << c << endl;
                        // cout << "       t :" << t1 << endl;
                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
            }

            condition = false;

            if (i == 0)
            {
                condition = !(spheres_list[i + 1]->isInside(step * t2 + O, limit_length));
            }
            else if (i == spheres_list.size() - 1)
            {
                condition = !(spheres_list[i - 1]->isInside(step * t2 + O, limit_length));
            }
            else
            {
                condition = !(spheres_list[i + 1]->isInside(step * t2 + O, limit_length) || spheres_list[i - 1]->isInside(step * t2 + O, limit_length));
                // cout <<  "spheres[" << i+1 << "].minDistance(step*t1+O) : "<<spheres[i+1].minDistance(step*t2+O) << endl;
                // cout <<  "spheres[" << i-1 << "].minDistance(step*t1+O) : "<<spheres[i-1].minDistance(step*t2+O) << endl;
                // cout << "this sphere : spheres[" << i << "].minDistance(step*t1+O) : "<<spheres[i].minDistance(step*t2+O) << endl;
                // cout << "condition : " << condition << endl;
            }
            // condition = true;
            // if the new position is at edge of axon, at the edge of sphere[i] but inside the neighbours
            if (condition)
            {
                if (Walker::bouncing)
                {
                    // You compare with esp_val to prevent taking the collision point into account for the next collision
                    if (t2 >= EPS_VAL && t2 <= step_lenght + barrier_tickness)
                    {
                        // cout << "   intersection, sphere :" << i << endl;
                        // cout << "       c :" << c<< endl;
                        // cout << "       t :" << t2 << endl;
                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
                else
                {
                    if (t2 >= 0 && t2 <= step_lenght + barrier_tickness)
                    {
                        // cout << "   intersection, sphere :" << i << endl;
                        // cout << "       c :" << c<< endl;
                        // cout << "       t :" << t2 << endl;
                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
            }
        }
    }

    if (dist_intersections.size() > 0)
    {
        // cout << "dist_intersections.size() :" << dist_intersections.size() << endl;

        auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
        unsigned index_ = std::distance(std::begin(dist_intersections), min_distance_int);

        int sphere_ind = sph_ids[index_];

        double dist_to_collision = dist_intersections[index_];

        colision.type = Collision::hit;
        colision.rn = cs[index_];

        if (walker.initial_location == Walker::intra)
        {
            if (check_inside_(all_cs))
            {
                colision.col_location = Collision::inside;
            }
            else
            {
                colision.col_location = Collision::outside;
            }
        }
        else if (walker.initial_location == Walker::extra)
        {
            if (check_outside_(all_cs))
            {
                colision.col_location = Collision::outside;
            }
            else
            {
                colision.col_location = Collision::inside;
            }
        }

        // cout << "Collision, sphere :" << sph_ids[index_] << endl;
        colision.t = fmin(dist_to_collision, step_lenght);
        cout << walker.pos_v[0] << " " << walker.pos_v[1] << " " << walker.pos_v[2] << endl;
        cout << "t " << colision.t << " step " << step[0] << " " << step[1] << " " << step[2] << endl;
        colision.colision_point = walker.pos_v + colision.t * step;
        cout << "t " << colision.colision_point[0] << " " << colision.colision_point[1] << " " << colision.colision_point[2] << endl;

        Vector3d normal = (colision.colision_point - spheres_list[sphere_ind]->center).normalized();

        Vector3d temp_step = step;
        elasticBounceAgainsPlane_intra(walker.pos_v, normal, colision.t, temp_step);
        colision.bounced_direction = temp_step.normalized();
        // colision.collision_objects = {sphere_ind, id};

        return true;
    }
    else
    {
        // cout << " ----- " << endl;
        // cout << "No collision" << endl;
        colision.type = Collision::null;
        return false;
    }
}

vector<int> Neuron::closest_subbranch(Vector3d const &position, int const &dendrite_id, int const &subbranch_id, double const &step_length)
{
    const auto &subbranch = dendrites[dendrite_id].subbranches[subbranch_id];
    int nb_subbranches = dendrites[dendrite_id].subbranches.size();
    int size_subbranch = subbranch.spheres.size() - 1;
    double eps = 0.0005;
    double distance_to_proximal_end = (subbranch.spheres[0].center - position).norm() - eps;
    double distance_to_distal_end = (subbranch.spheres[size_subbranch].center - position).norm() - eps;

    // cout << "d prox" << distance_to_proximal_end << endl;
    // cout << "d dist" << distance_to_distal_end << endl;
    if (distance_to_proximal_end <= step_length)
    {
        // The subbranch id starts at 1 but the indices in c++ start at 0
        vector<int> proximal_branching = subbranch.proximal_branching;
        if (proximal_branching.size() > 1)
        {
            proximal_branching[0] -= 1;
            proximal_branching[1] -= 1;
            //    cout << "p" << proximal_branching[0] << " " << proximal_branching[1] << endl;
        }
        return proximal_branching;
    }
    else if (distance_to_distal_end <= step_length)
    {
        // The subbranch id starts at 1 but the indices in c++ start at 0
        vector<int> distal_branching = subbranch.distal_branching;

        // No need to check the size, because there are always 2 distal branching
        distal_branching[0] -= 1;
        distal_branching[1] -= 1;
        // cout << "d" << distal_branching[0] << " " << distal_branching[1] << endl;

        // Distal end of the neuron
        if ((distal_branching[0] >= nb_subbranches) || (distal_branching[1] >= nb_subbranches))
            return {-2};

        return distal_branching;
    }
    else
        return {-1};
}

int Neuron::closest_dendrite_from_soma(Vector3d const &position, double const &step_length)
{
    int smaller_dist = 1000;
    int closer_dendrite = -1;
    for (size_t i = 0; i < dendrites.size(); ++i)
    {
        // TODO : there was a segfault here once ... [ines]
        double distance_to_axon = (position - dendrites[i].subbranches[0].spheres[0].center).norm();
        if (distance_to_axon < smaller_dist && distance_to_axon <= step_length)
        {
            smaller_dist = distance_to_axon;
            closer_dendrite = i;
        }
    }
    return closer_dendrite;
}

// TODO: do other cases
TEST_CASE("intersection_sphere_vector")
{
    Vector3d center(0.05, 0.05, 0.05);
    double radius = 5e-4;
    int sphere_id = 0;
    Dynamic_Sphere sphere(center, radius, sphere_id);

    Vector3d direction(1, 0, 0);
    double step_length = 5e-4;

    Vector3d traj_origin(0.05, 0.05, 0.05);

    double intercept1, intercept2, c;
    CHECK(Neuron::intersection_sphere_vector(intercept1, intercept2, sphere, direction, step_length, traj_origin, c));

    CHECK_EQ(intercept1, doctest::Approx(radius));
    CHECK_EQ(intercept2, doctest::Approx(-radius));
    // CHECK_EQ(c, ); // TODO: [ines]
}

bool Neuron::intersection_sphere_vector(double &intercept1, double &intercept2, Dynamic_Sphere const &sphere, Vector3d const &step_dir, double const &step_length, Vector3d const &traj_origin, double &c)
{
    // cout << sphere.center[0] << endl;
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

void Neuron::set_dendrites(std::vector<Dendrite> const &dendrites_to_add)
{
    if (dendrites_to_add.size() != 0)
        dendrites = dendrites_to_add;
}

void Neuron::add_dendrite(Dendrite &dendrite_to_add)
{
    dendrite_to_add.add_projection(dendrite_to_add.id);
    dendrites.push_back(dendrite_to_add);
}

void Neuron::generateSpanRadius(double const &lower_bound, double const &upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<double> dist_span_radius(lower_bound, upper_bound);
    // cout << dist_span_radius(rng) << endl;
    // Generate int number in [lb, ub], in [mm]
    span_radius = dist_span_radius(rng);
}

vector<double> Neuron::get_Volume() const
{
    double VolumeSoma = 0;
    double VolumeDendrites = 0;
    // Calculate the volume of the soma
    VolumeSoma += 4 / 3 * M_PI * pow(soma.radius, 3);

    // Calculate the cylindrical volume of each dendrite
    for (uint8_t j = 0; j < dendrites.size(); j++)
    {
        VolumeDendrites += dendrites[j].volumeDendrite();
    }

    return {VolumeSoma, VolumeDendrites};
}