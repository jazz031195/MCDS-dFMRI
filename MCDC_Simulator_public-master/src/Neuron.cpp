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
    nb_dendrites = 20; // dist_dendrites(rng);

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

TEST_CASE("isPosInsideNeuron")
{
    cout << "isPosInsideNeuron" << endl;
    Vector3d center(0.05, 0.05, 0.05);
    double radius_soma     = 10e-3;
    double radius_dendrite = 0.5e-3;
    Neuron neuron(center, radius_soma, 0);
    Dendrite dendrite;

    int branch_id = 0;
    vector<Dynamic_Sphere> spheres_list;
    for (size_t i = 0; i < 10; ++i)
    {
        Vector3d next_center(center[0] + radius_soma + i * radius_dendrite / 4, center[1], center[2]);
        Dynamic_Sphere sphere_to_add(next_center, radius_dendrite, 0, false, branch_id, i, 1);
        spheres_list.push_back(sphere_to_add);
    }

    Vector3d last_center(spheres_list[spheres_list.size()-1].center);
    Vector3d begin;
    vector <int> proximal_end = {0};
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
        branching_direction = {1, -1, 0};
        branching_direction = branching_direction.normalized();
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
    // vector<int> sphere_ids;
    // for (size_t p = 0; p < dendrites.size(); p++)
    // {
    //     for (size_t b = 0; b < dendrites[p].subbranches.size(); b++)
    //     {
    //         for (size_t s = 0; s < dendrites[p].subbranches[b].spheres.size(); s++)
    //         {
    //             Vector3d somaCenter = dendrites[p].subbranches[b].spheres[s].center;
    //             double somaRadius = dendrites[p].subbranches[b].spheres[s].radius;

    //             if (pow(position[0] - somaCenter[0], 2) + pow(position[1] - somaCenter[1], 2) + pow(position[2] - somaCenter[2], 2) <= somaRadius + distance_to_be_inside)
    //             {
    //                 in_dendrite_index = p;
    //                 in_subbranch_index = b;
    //                 in_soma_index = -1;
    //                 in_sph_index = sphere_ids;
    //                 return true;
    //             }
    //         }
    //     }
    // }

    // if (soma.isInside(position, distance_to_be_inside))
    // {
    //     in_soma_index = 0;
    //     in_dendrite_index = -1;
    //     in_subbranch_index = -1;
    //     in_sph_index.clear();
    //     return true;
    // }

    // if position is in box with Dendrite inside
    // distance_to_be_inside = position.sphere.radius + constant
    // id of the dendrite. {} if not in neuron
    vector<int> part_id = isNearDendrite(position, distance_to_be_inside);
    // // cout << neuron_part << endl;
    vector<int> sphere_ids;
    if (part_id.size() > 0)
    {
        for (size_t p=0; p < part_id.size(); p++)
        {
            for (size_t b=0; b < dendrites[part_id[p]].subbranches.size(); b++)
            {
                if (dendrites[part_id[p]].subbranches[b].isPosInsideAxon(position, distance_to_be_inside, sphere_ids))
                {
                    in_dendrite_index  = part_id[p];
                    in_subbranch_index = b;
                    in_soma_index      = -1;
                    in_sph_index       = sphere_ids;
                    return true;
                }
            }
        }
    }

    if (isNearSoma(position, distance_to_be_inside))
    {
        if (soma.isInside(position, distance_to_be_inside))
        {
            in_soma_index      = 0;
            in_dendrite_index  = -1;
            in_subbranch_index = -1;
            in_sph_index.clear();
            return true;
        }
    }
    in_dendrite_index  = -1;
    in_subbranch_index = -1;
    in_soma_index      = -1;
    in_sph_index.clear();
    return false;
}

TEST_CASE("isNearSoma")
{
    cout << "isNearSoma" << endl;
    Vector3d center(0.05, 0.05, 0.05);
    double radius = 10e-3;
    int sphere_id = 0;
    Neuron neuron(center, radius, sphere_id);

    SUBCASE("inside soma")
    {
        Vector3d position(center[0] + radius - 2*EPS_VAL, center[1], center[2]);
        CHECK(neuron.isNearSoma(position, EPS_VAL));
    }
    SUBCASE("outside soma")
    {
        Vector3d position(center[0] + radius + 2*EPS_VAL, center[1], center[2]);
        CHECK(!neuron.isNearSoma(position, EPS_VAL));
    }
}

bool Neuron::isNearSoma(Vector3d const &position, double const &distance_to_be_inside) const
{
    int count_isnear = 0;
    // Check soma box
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
        if ((position[axis] >= soma.center[axis] - soma.radius - distance_to_be_inside) &&
            (position[axis] <= soma.center[axis] + soma.radius + distance_to_be_inside))
        {
            ++count_isnear;
        }
    }
    if (count_isnear == 3)
        return true;
    else
        return false;
}

TEST_CASE("isNearDendrite")
{
    cout << "isNearDendrite" << endl;
    Vector3d center(0.05, 0.05, 0.05);
    double radius_soma = 10e-3;
    double radius_dendrite = 0.5e-3;

    Neuron neuron(center, radius_soma, 0);
    Dendrite dendrite;

    int branch_id = 0;
    vector<Dynamic_Sphere> spheres_list;
    for (size_t i = 0; i < 11; ++i)
    {
        Vector3d next_center(center[0] + radius_soma + i * radius_dendrite / 4, center[1], center[2]);
        Dynamic_Sphere sphere_to_add(next_center, radius_dendrite, 0, false, branch_id, i, 1);
        spheres_list.push_back(sphere_to_add);
    }

    Vector3d begin;
    vector <int> proximal_end = {1};
    vector <int> distal_end   = {2, 3};
    Axon subbranch(branch_id, radius_dendrite, begin, begin, 0, false, 1, proximal_end, distal_end);
    subbranch.set_spheres(spheres_list);
    dendrite.add_subbranch(subbranch);
    neuron.add_dendrite(dendrite);

    spheres_list.clear();
    for (size_t i = 0; i < 11; ++i)
    {
        Vector3d next_center(center[0] - radius_soma + i * radius_dendrite / 4, center[1], center[2]);
        Dynamic_Sphere sphere_to_add(next_center, radius_dendrite, 0, false, branch_id, i, 1);
        spheres_list.push_back(sphere_to_add);
    }

    Dendrite dendrite2;
    Axon subbranch2(branch_id + 1, radius_dendrite, begin, begin, 0, false, 1, proximal_end, distal_end);
    subbranch2.set_spheres(spheres_list);
    dendrite2.add_subbranch(subbranch2);
    neuron.add_dendrite(dendrite2);

    vector<int> dendrite_ids;
    SUBCASE("Soma center")
    {
        Vector3d position(center[0], center[1], center[2]);
        dendrite_ids = neuron.isNearDendrite(position, EPS_VAL);
        // Should not be near any dendrite
        CHECK(dendrite_ids.size() == 0);
    }
    SUBCASE("Begin of dendrite 0")
    {
        Vector3d position(center[0] + radius_soma + 2*EPS_VAL, center[1], center[2]);
        dendrite_ids = neuron.isNearDendrite(position, EPS_VAL);
        CHECK(dendrite_ids.size() == 1);
        CHECK_EQ(dendrite_ids[0], 0);
    }
    SUBCASE("End of dendrite 0")
    {
        auto sphere_list = neuron.dendrites[0].subbranches[0].spheres;
        Vector3d position(sphere_list[sphere_list.size()-1].center);
        dendrite_ids = neuron.isNearDendrite(position, EPS_VAL);
        CHECK(dendrite_ids.size() == 1);
        CHECK_EQ(dendrite_ids[0], 0);
    }
    SUBCASE("Begin: dendrite 1")
    {
        Vector3d position(center[0] - radius_soma - 2*EPS_VAL, center[1], center[2]);
        dendrite_ids = neuron.isNearDendrite(position, EPS_VAL);
        CHECK(dendrite_ids.size() == 1);
        CHECK_EQ(dendrite_ids[0], 1);
    }
    SUBCASE("End: dendrite 1")
    {
        auto sphere_list = neuron.dendrites[1].subbranches[0].spheres;
        Vector3d position(sphere_list[sphere_list.size()-1].center);
        dendrite_ids = neuron.isNearDendrite(position, EPS_VAL);
        CHECK(dendrite_ids.size() == 1);
        CHECK_EQ(dendrite_ids[0], 1);
    }
}

vector<int> Neuron::isNearDendrite(Vector3d const &position, double const &distance_to_be_inside) const
{
    // Check each dendrite's box
    int count_isnear = 0;
    vector<int> dendrite_ids;
    for (unsigned int i = 0; i < dendrites.size(); ++i)
    {
        count_isnear = 0;
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
            Vector2d axis_limits = dendrites[i].projections.axon_projections[axis];
            if ((position[axis] >= axis_limits[0] - distance_to_be_inside) &&
                (position[axis] <= axis_limits[1] + distance_to_be_inside))
            {
                ++count_isnear;
            }
        }
        // Inside the box around dendrite
        if (count_isnear == 3)
            dendrite_ids.push_back(i);
    }
    if (dendrite_ids.size() > 0)
        return dendrite_ids;

    return dendrite_ids;
}

TEST_CASE("find_neighboring_spheres")
{
    cout << "find_neighboring_spheres" << endl;
    Vector3d center(0.05, 0.05, 0.05);
    double radius_soma     = 10e-3;
    double radius_dendrite = 0.5e-3;
    Neuron neuron(center, radius_soma, 0);
    Dendrite dendrite;

    int branch_id = 0;
    vector<Dynamic_Sphere> spheres_list;
    for (size_t i = 0; i < 20; ++i)
    {
        Vector3d next_center(center[0] + radius_soma + i * radius_dendrite / 4, center[1], center[2]);
        Dynamic_Sphere sphere_to_add(next_center, radius_dendrite, 0, false, branch_id, i, 1);
        spheres_list.push_back(sphere_to_add);
    }

    Vector3d last_center(spheres_list[spheres_list.size()-1].center);
    Vector3d begin;
    vector <int> proximal_end = {0};
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
        for (size_t i = 1; i < 20; ++i)
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
        branching_direction = {1, -1, 0};
        branching_direction = branching_direction.normalized();
    }
    neuron.add_dendrite(dendrite);

    Walker walker;
    auto empty_vect = vector<unsigned int>();
    walker.dyn_cylinders_collision_sphere.collision_list = &empty_vect;
    walker.cylinders_collision_sphere.collision_list     = &empty_vect;
    walker.spheres_collision_sphere.collision_list       = &empty_vect;
    walker.dyn_spheres_collision_sphere.collision_list   = &empty_vect;
    walker.axons_collision_sphere.collision_list         = &empty_vect;
    auto empty_vect_vect = vector<vector<unsigned int>>();
    walker.ply_collision_sphere.collision_list           = &empty_vect_vect;

    auto neur_list = vector<unsigned int>({0});
    walker.neurons_collision_sphere.collision_list = &neur_list;
    walker.neurons_collision_sphere.pushToBigSphere(0);
    walker.neurons_collision_sphere.pushToSmallSphere(0);
    walker.in_neuron_index = 0;

    double step_length = 5.5e-4; 

    SUBCASE("Soma center")
    {
        auto ini_pos = Vector3d(center[0], center[1], center[2]);
        walker.setInitialPosition(ini_pos);
        Vector3d step_dir(1, 0, 0);
        neuron.isPosInsideNeuron(ini_pos, barrier_tickness, false, walker.in_soma_index, walker.in_dendrite_index,
                                walker.in_subbranch_index, walker.in_sph_index);
        vector<vector<Dynamic_Sphere *>> spheres_list = neuron.find_neighbor_spheres(walker, step_dir, step_length);
        CHECK_EQ(spheres_list.size(), 1);
        CHECK_EQ(spheres_list[0].size(), 1);
        CHECK_EQ(spheres_list[0][0]->radius, 10e-3);
    }
    SUBCASE("Soma radius")
    {
        auto ini_pos = Vector3d(center[0] + radius_soma - 2*radius_dendrite, center[1], center[2]);
        walker.setInitialPosition(ini_pos);
        Vector3d step_dir(1, 0, 0);
        neuron.isPosInsideNeuron(ini_pos, barrier_tickness, false, walker.in_soma_index, walker.in_dendrite_index,
                                walker.in_subbranch_index, walker.in_sph_index);
        vector<vector<Dynamic_Sphere *>> spheres_list = neuron.find_neighbor_spheres(walker, step_dir, step_length);
        CHECK_EQ(spheres_list.size(), 2);
        CHECK_EQ(spheres_list[0].size(), 1);
        CHECK_EQ(spheres_list[1].size(), int(step_length/radius_dendrite*4));
        CHECK_EQ(spheres_list[0][0]->radius, 10e-3);
        CHECK_EQ(spheres_list[0][0]->id, 0);
        CHECK_EQ(spheres_list[1][0]->radius, 0.5e-3);
        CHECK_EQ(spheres_list[1][0]->id, 0);
    }
    SUBCASE("Branching")
    {
        SUBCASE("End subbranch 0")
        {
            auto ini_pos = neuron.dendrites[0].subbranches[0].spheres[19].center;
            walker.setInitialPosition(ini_pos);
            Vector3d step_dir(1, 0, 0);
            neuron.isPosInsideNeuron(ini_pos, barrier_tickness, false, walker.in_soma_index, walker.in_dendrite_index,
                                    walker.in_subbranch_index, walker.in_sph_index);
            vector<vector<Dynamic_Sphere *>> spheres_list = neuron.find_neighbor_spheres(walker, step_dir, step_length);
            CHECK_EQ(spheres_list.size(), 3);
            CHECK_EQ(spheres_list[0][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[1][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[1][0]->id, 0);
            CHECK_EQ(spheres_list[1].size(), 10);
            CHECK_EQ(spheres_list[2][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[2][0]->id, 0);
            CHECK_EQ(spheres_list[2].size(), 10);
        }
        SUBCASE("Begin subbranch 1")
        {
            auto ini_pos = neuron.dendrites[0].subbranches[1].spheres[0].center;
            walker.setInitialPosition(ini_pos);
            Vector3d step_dir(1, -1, 0);
            step_dir = step_dir.normalized();
            neuron.isPosInsideNeuron(ini_pos, barrier_tickness, false, walker.in_soma_index, walker.in_dendrite_index,
                                    walker.in_subbranch_index, walker.in_sph_index);
            vector<vector<Dynamic_Sphere *>> spheres_list = neuron.find_neighbor_spheres(walker, step_dir, step_length);
            CHECK_EQ(spheres_list.size(), 3);
            CHECK_EQ(spheres_list[0][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[1][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[1][0]->id, 0);
            CHECK_EQ(spheres_list[1].size(), 10);
            CHECK_EQ(spheres_list[2][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[2][0]->id, 0);
            CHECK_EQ(spheres_list[2].size(), 10);
        }
        SUBCASE("Begin subbranch 2")
        {
            auto ini_pos = neuron.dendrites[0].subbranches[2].spheres[0].center;
            walker.setInitialPosition(ini_pos);
            Vector3d step_dir(-1, 1, 0);
            step_dir = step_dir.normalized();
            neuron.isPosInsideNeuron(ini_pos, barrier_tickness, false, walker.in_soma_index, walker.in_dendrite_index,
                                    walker.in_subbranch_index, walker.in_sph_index);
            vector<vector<Dynamic_Sphere *>> spheres_list = neuron.find_neighbor_spheres(walker, step_dir, step_length);
            CHECK_EQ(spheres_list.size(), 3);
            CHECK_EQ(spheres_list[0][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[1][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[1][0]->id, 0);
            CHECK_EQ(spheres_list[1].size(), 10);
            CHECK_EQ(spheres_list[2][0]->radius, 0.5e-3);
            CHECK_EQ(spheres_list[2][0]->id, 0);
            CHECK_EQ(spheres_list[2].size(), 10);
            CHECK_EQ((spheres_list[0][spheres_list[0].size()-1]->center - spheres_list[1][0]->center).norm(),  doctest::Approx(radius_dendrite/4));
            CHECK_EQ((spheres_list[0][spheres_list[0].size()-1]->center - spheres_list[2][0]->center).norm(),  doctest::Approx(radius_dendrite/4));
        }
    
    }

}

vector<vector<Dynamic_Sphere *>> Neuron::find_neighbor_spheres(Walker &walker, Vector3d const &step_dir, double const &step_length)
{
    vector<vector<Dynamic_Sphere *>> spheres_list;
    vector<Dynamic_Sphere *> spheres_sublist1, spheres_sublist2, spheres_sublist3;

    // spheres_list.push_back(new Dynamic_Sphere(soma));

    // if (walker.in_soma_index == 0)
    // {
    //     return spheres_list;
    // }
    // if ((walker.in_dendrite_index == -1) || (walker.in_subbranch_index == -1))
    // {
    //     cout << "weird" << endl;
    //     return spheres_list;
    // }

    // const auto subbranch = dendrites[walker.in_dendrite_index].subbranches[walker.in_subbranch_index];
    // for (size_t i = 0; i < subbranch.spheres.size(); i++)
    // {
    //     spheres_list.push_back(new Dynamic_Sphere(subbranch.spheres[i]));
    // }

    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d next_pt = O + step_length * step_dir;
    // We are in the soma
    if(walker.in_soma_index == 0)
    {
        spheres_sublist1.push_back(new Dynamic_Sphere(soma));
       // Find in which dendrite it can go
        int closest_dendrite_id = closest_dendrite_from_soma(next_pt, step_length);

        // If inside closest_dendrite_id dendrite => no collision
        // TODO : check if enough to check only the first subbranch [ines]
        if (closest_dendrite_id >= 0)
        {
            size_t nb_spheres_needed = size_t(step_length / (dendrites[closest_dendrite_id].subbranches[0].spheres[0].radius/4));
            for(size_t i=0; i < nb_spheres_needed; ++i)
            {
                spheres_sublist2.push_back(new Dynamic_Sphere(dendrites[closest_dendrite_id].subbranches[0].spheres[i]));
            }
        }
        spheres_list.push_back(spheres_sublist1);
        if(spheres_sublist2.size() > 0)
            spheres_list.push_back(spheres_sublist2);
    }
    else if(walker.in_dendrite_index >= 0)
    {
        // Calculate the dot product between step_dir and (subbranch.spheres[end-1].center - O)
        auto current_subbranch   = dendrites[walker.in_dendrite_index].subbranches[walker.in_subbranch_index];
        size_t current_sphere_first = size_t(walker.in_sph_index[0]);
        size_t current_sphere_last  = size_t(walker.in_sph_index[walker.in_sph_index.size() - 1]);

        // Goes to proximal
        if(current_sphere_first < 10)
        {
            //check proximal branching
            vector<int> proximal_branching = current_subbranch.proximal_branching;
            // Close to the soma
            if(proximal_branching.size() == 1)
                spheres_sublist1.push_back(new Dynamic_Sphere(soma));
            else
            {
                for(size_t proximal_id=0; proximal_id < 10; ++ proximal_id)
                {
                    auto proximal_branch1_spheres = dendrites[walker.in_dendrite_index].subbranches[proximal_branching[0] - 1].spheres;
                    if(proximal_branch1_spheres.size() > proximal_id)
                        spheres_sublist1.push_back(new Dynamic_Sphere(proximal_branch1_spheres[proximal_branch1_spheres.size() - 10 + proximal_id]));
                    auto proximal_branch2_spheres = dendrites[walker.in_dendrite_index].subbranches[proximal_branching[1] - 1].spheres;
                    if(proximal_branch2_spheres.size() > proximal_id)
                        spheres_sublist3.push_back(new Dynamic_Sphere(proximal_branch2_spheres[proximal_id]));
                }
            }
            // Add the sphere in which the walker is
            for(size_t sph_id=0; sph_id < walker.in_sph_index.size() + 10; ++ sph_id)
                spheres_sublist2.push_back(new Dynamic_Sphere(current_subbranch.spheres[sph_id]));

            if(spheres_sublist1.size() > 0)
                spheres_list.push_back(spheres_sublist1);
            if(spheres_sublist2.size() > 0)
                spheres_list.push_back(spheres_sublist2);
            if(spheres_sublist3.size() > 0)
                spheres_list.push_back(spheres_sublist3);
        }
        // Goes to distal
        else if(current_sphere_last + 10 > current_subbranch.spheres.size())
        {
            // Add the sphere in which the walker is
            for(size_t sph_id=current_sphere_first - 10; sph_id < current_subbranch.spheres.size(); ++ sph_id)
                spheres_sublist1.push_back(new Dynamic_Sphere(current_subbranch.spheres[sph_id]));

            //check distal branching
            vector<int> distal_branching = current_subbranch.distal_branching;
            if(distal_branching[0] < dendrites[walker.in_dendrite_index].subbranches.size())
            {
                for(size_t distal_id=0; distal_id < 10; ++ distal_id)
                {
                    auto distal_branch1_spheres = dendrites[walker.in_dendrite_index].subbranches[distal_branching[0] - 1].spheres;
                    if(distal_branch1_spheres.size() > distal_id)
                        spheres_sublist2.push_back(new Dynamic_Sphere(distal_branch1_spheres[distal_id]));
                    auto proximal_branch2_spheres = dendrites[walker.in_dendrite_index].subbranches[distal_branching[1] - 1].spheres;
                    if(proximal_branch2_spheres.size() > distal_id)
                        spheres_sublist3.push_back(new Dynamic_Sphere(proximal_branch2_spheres[distal_id]));
                }
            }
            if(spheres_sublist1.size() > 0)
                spheres_list.push_back(spheres_sublist1);
            if(spheres_sublist2.size() > 0)
                spheres_list.push_back(spheres_sublist2);
            if(spheres_sublist3.size() > 0)
                spheres_list.push_back(spheres_sublist3); 
        }
        else
        {
            // Add the sphere in which the walker is
            for(size_t sph_id=current_sphere_first-10; sph_id < current_sphere_last + 10; ++ sph_id)
                spheres_sublist1.push_back(new Dynamic_Sphere(current_subbranch.spheres[sph_id]));

            spheres_list.push_back(spheres_sublist1);
        }               
    }
    return spheres_list;
}

bool Neuron::checkCollision(Walker &walker, Vector3d const &step_dir, double const &step_lenght, Collision &colision)
{
    bool isColliding = false;

    vector<vector<Dynamic_Sphere *>> spheres_list = find_neighbor_spheres(walker, step_dir, step_lenght);

    isColliding = checkCollision_branching(walker, spheres_list, step_dir, step_lenght, colision);


    for (size_t i = 0; i < spheres_list.size(); i++)
        for(size_t j = 0; j < spheres_list[i].size(); j++)
            delete spheres_list[i][j];
    spheres_list.clear();

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
    cout << "checkCollision_branching" << endl;
    Vector3d center(0.05, 0.05, 0.05);
    double radius_soma = 10e-3;
    double radius_dendrite = 0.5e-3;
    int sphere_id = 0;

    vector<vector<Dynamic_Sphere *>> spheres_list;
    
    SUBCASE("Soma + 1 branch")
    {
        spheres_list.push_back({new Dynamic_Sphere(center, radius_soma, sphere_id)});
        vector<Dynamic_Sphere *> spheres_list_tmp;
        for (size_t i = 0; i < 20; ++i)
        {
            Vector3d next_center(center[0] + radius_soma + i * radius_dendrite / 4, center[1], center[2]);
            spheres_list_tmp.push_back(new Dynamic_Sphere(next_center, radius_dendrite, i));
        }
        spheres_list.push_back(spheres_list_tmp);
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
            w.in_soma_index = 0;

            bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
            CHECK(!collided);
            CHECK_EQ(colision.type, Collision::null);
            CHECK_EQ(colision.t, INFINITY_VALUE);
        }
        SUBCASE("no collision: soma to dendrite")
        {
            Vector3d direction(1, 0, 0);
            double step_length = 5e-4;

            Walker w;
            // w.setNextDirection(direction);
            w.setInitialPosition(center[0] + radius_soma - step_length/2, center[1], center[2]);
            w.location = w.initial_location = Walker::intra;
            w.in_soma_index = 0;

            bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
            CHECK(!collided);
            CHECK_EQ(colision.type, Collision::null);
            CHECK_EQ(colision.t, INFINITY_VALUE);
        }
        SUBCASE("no collision: dendrite to soma")
        {
            Vector3d direction(-1, 0, 0);
            double step_length = 5e-4;

            Walker w;
            Vector3d center_first_dendrite(spheres_list[1][0]->center);
            w.setInitialPosition(center_first_dendrite[0] + step_length/2, center_first_dendrite[1], center_first_dendrite[2]);
            w.location = w.initial_location = Walker::intra;
            w.in_soma_index = 0;

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
            w.in_soma_index = 0;

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
    }
    SUBCASE("Branching")
    {
        int branch_id = 0;
        vector<Dynamic_Sphere*> spheres_list_tmp1;
        for (size_t i = 0; i < 10; ++i)
        {
            Vector3d next_center(center[0] + i * radius_dendrite / 4, center[1], center[2]);
            spheres_list_tmp1.push_back(new Dynamic_Sphere(next_center, radius_dendrite, 0, false, branch_id, i, 1));
        }
        spheres_list.push_back(spheres_list_tmp1);
        Vector3d last_center(spheres_list_tmp1[spheres_list_tmp1.size()-1]->center);

        // Add branching
        Vector3d branching_direction(1, 1, 0);
        branching_direction = branching_direction.normalized();

        for(size_t rep=0; rep < 2; rep++)
        {
            spheres_list_tmp1.clear();
            for (size_t i = 1; i < 11; ++i)
            {
                Vector3d next_center(last_center[0] + i * radius_dendrite / 4 * branching_direction[0],
                                     last_center[1] + i * radius_dendrite / 4 * branching_direction[1],
                                     last_center[2] + i * radius_dendrite / 4 * branching_direction[2]);
                spheres_list_tmp1.push_back(new Dynamic_Sphere(next_center, radius_dendrite, 0, false, branch_id + rep + 1, i, 1));
            }
            spheres_list.push_back(spheres_list_tmp1);
            branching_direction = {1, -1, 1};
        }

        Collision colision;
        colision.type = Collision::null;
        colision.t = INFINITY_VALUE;

        double step_length = 7e-4;

        Walker w;
        SUBCASE("proximal to distal")
        {   
            Vector3d center = spheres_list[0][spheres_list[0].size() - 1]->center;
            w.setInitialPosition(center[0], center[1], center[2]);
            w.location = w.initial_location = Walker::intra;
            w.in_dendrite_index  = 0;
            w.in_subbranch_index = 0;
            w.in_sph_index       = {int(spheres_list[0].size() - 3), int(spheres_list[0].size() - 2), int(spheres_list[0].size() - 1)};

            SUBCASE("proximal to distal 1")
            { 
                Vector3d direction(1, 1, 0);
                direction = direction.normalized();

                bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
                CHECK(!collided);
                CHECK_EQ(colision.type, Collision::null);
            }
            SUBCASE("proximal to distal 2")
            {   
                Vector3d direction(1, -1, 0);
                direction = direction.normalized();

                bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
                CHECK(!collided);
                CHECK_EQ(colision.type, Collision::null);
            }
        }
        SUBCASE("Distal to proximal 1")
        {
            Vector3d center = spheres_list[1][0]->center;
            w.setInitialPosition(center[0], center[1], center[2]);
            w.location = w.initial_location = Walker::intra;
            w.in_dendrite_index  = 0;
            w.in_subbranch_index = 1;
            w.in_sph_index       = {0, 1, 2};

            Vector3d direction(-1, -1, 0);
            direction = direction.normalized();

            bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
            CHECK(!collided);
            CHECK_EQ(colision.type, Collision::null);
        }
        SUBCASE("Distal to proximal 2")
        {
            Vector3d center = spheres_list[2][0]->center;
            w.setInitialPosition(center[0], center[1], center[2]);
            w.location = w.initial_location = Walker::intra;
            w.in_dendrite_index  = 0;
            w.in_subbranch_index = 1;
            w.in_sph_index       = {0, 1, 2};

            Vector3d direction(-1, 1, 0);
            direction = direction.normalized();

            bool collided = Neuron::checkCollision_branching(w, spheres_list, direction, step_length, colision);
            CHECK(!collided);
            CHECK_EQ(colision.type, Collision::null);
        }
    
    }
    
    // TODO: [ines] more complex cases
    for (size_t i = 0; i < spheres_list.size(); i++)
        for(size_t j = 0; j < spheres_list[i].size(); j++)
            delete spheres_list[i][j];
    spheres_list.clear();
}

bool Neuron::checkCollision_branching(Walker &walker, vector<vector<Dynamic_Sphere *>> &spheres_list, Vector3d const &step, double const &step_lenght, Collision &colision)
{

    string message;
    Vector3d O;
    walker.getVoxelPosition(O);

    // distances to intersections
    std::vector<double> dist_intersections;
    // values indicating whether the walker is inside or outside a sphere
    std::vector<double> cs;

    std::vector<int> sph_ids;
    std::vector<int> branch_ids;
    std::vector<double> all_cs;

    for(size_t list_id=0; list_id < spheres_list.size(); ++list_id)
    {
        for (size_t i = 0; i < spheres_list[list_id].size(); ++i)
        {

            // distances to collision
            double t1;
            double t2;
            double c;
            bool intersect = intersection_sphere_vector(t1, t2, *spheres_list[list_id][i], step, step_lenght, O, c);
            double limit_length = -EPS_VAL;

            if (intersect)
            {
                all_cs.push_back(c);
                
                // if the new position is at edge of axon, at the edge of sphere[i] but inside the neighbours
                bool condition = false;

                // Checking the first list
                if(list_id == 0)
                {
                    // Colliding with soma only
                    if((spheres_list[list_id].size() == 1) && (spheres_list.size() == 1))
                    {
                        condition = true;
                    }
                    // In the first branch
                    else if ((i == 0) && (spheres_list[list_id].size() > 1))
                    {
                        condition = !(spheres_list[list_id][i + 1]->isInside(step * t1 + O, limit_length));
                    }
                    // We are at the distal end of the branch
                    else if(i == spheres_list[list_id].size() - 1)
                    {
                        // There is a sphere before
                        if (spheres_list[list_id].size() > 1)
                        {
                            condition = (spheres_list[list_id][i - 1]->isInside(step * t1 + O, limit_length));
                        }
                        // There is a distal branch 
                        if(spheres_list.size() > 1)
                        {
                            for(size_t dist_branch=1; dist_branch < spheres_list.size(); dist_branch++)
                            {
                                condition = condition || (spheres_list[dist_branch][0]->isInside(step * t1 + O, limit_length));
                            }
                        }
                        condition = !condition;
                    }
                    // We are in the middle of the branch
                    else
                    {   
                        condition = !(spheres_list[list_id][i + 1]->isInside(step * t1 + O, limit_length) || spheres_list[list_id][i - 1]->isInside(step * t1 + O, limit_length));
                    }

                    if (condition)
                    {

                        // if the collision are too close or negative.
                        if (Walker::bouncing)
                        {
                            if (t1 >= EPS_VAL && t1 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t1);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                        else
                        {
                            if (t1 >= 0 && t1 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t1);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                    }

                    condition = false;

                    // Colliding with soma only
                    if((spheres_list[list_id].size() == 1) && (spheres_list.size() == 1))
                    {
                        condition = true;
                    }
                    // In the first branch
                    else if ((i == 0) && (spheres_list[list_id].size() > 1))
                    {
                        condition = !(spheres_list[list_id][i + 1]->isInside(step * t2 + O, limit_length));
                    }
                    // We are at the distal end of the branch
                    else if(i == spheres_list[list_id].size() - 1)
                    {
                        // There is a sphere before
                        if (spheres_list[list_id].size() > 1)
                        {
                            condition = (spheres_list[list_id][i - 1]->isInside(step * t2 + O, limit_length));
                        }
                        // There is a distal branch 
                        if(spheres_list.size() > 1)
                        {
                            for(size_t dist_branch=1; dist_branch < spheres_list.size(); dist_branch++)
                                condition = condition || (spheres_list[dist_branch][0]->isInside(step * t2 + O, limit_length));
                        }
                        condition = !condition;
                    }
                    // We are in the middle of the branch
                    else
                    {
                        condition = !(spheres_list[list_id][i + 1]->isInside(step * t2 + O, limit_length) || spheres_list[list_id][i - 1]->isInside(step * t2 + O, limit_length));
                    }

                    if (condition)
                    {

                        // if the collision are too close or negative.
                        if (Walker::bouncing)
                        {
                            if (t2 >= EPS_VAL && t2 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t2);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                        else
                        {
                            if (t2 >= 0 && t2 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t2);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                    }
                }
                // we are not in the first branch
                else
                {
                    // We should check the proximal branches
                    if(i == 0)
                    {
                        // Check the last sphere of the first proximal branch
                        condition = (spheres_list[0][spheres_list[0].size() - 1]->isInside(step * t1 + O, limit_length));

                        // There is another proximal branch
                        if(spheres_list.size() > 2)
                        {
                            condition = condition || (spheres_list[3-list_id][0]->isInside(step * t1 + O, limit_length));
                        }

                        // Check the sphere after
                        condition = !(condition || (spheres_list[list_id][i + 1]->isInside(step * t1 + O, limit_length)));

                    }
                    // We are at the end of the branch => check only the sphere before
                    else if(i == spheres_list[list_id].size() - 1)
                        condition = !(spheres_list[list_id][i - 1]->isInside(step * t1 + O, limit_length));
                    else
                    {
                        condition = !(spheres_list[list_id][i - 1]->isInside(step * t1 + O, limit_length) || spheres_list[list_id][i + 1]->isInside(step * t1 + O, limit_length));
                    }

                    if (condition)
                    {

                        // if the collision are too close or negative.
                        if (Walker::bouncing)
                        {
                            if (t1 >= EPS_VAL && t1 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t1);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                        else
                        {
                            if (t1 >= 0 && t1 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t1);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                    }
                    
                    condition = false;

                    // We should check the proximal branches
                    if(i == 0)
                    {
                        // Check the last sphere of the first proximal branch
                        condition = !(spheres_list[0][spheres_list[0].size() - 1]->isInside(step * t2 + O, limit_length));

                        // There is another proximal branch
                        if(spheres_list.size() > 2)
                        {
                            condition = condition || (spheres_list[3-list_id][0]->isInside(step * t2 + O, limit_length));
                        }

                        // Check the sphere after
                        condition = !(condition || (spheres_list[list_id][i + 1]->isInside(step * t2 + O, limit_length)));
                    }
                    // We are at the end of the branch => check only the sphere before
                    else if(i == spheres_list[list_id].size() - 1)
                        condition = !(spheres_list[list_id][i - 1]->isInside(step * t2 + O, limit_length));
                    else
                    {
                        condition = !(spheres_list[list_id][i - 1]->isInside(step * t2 + O, limit_length) || spheres_list[list_id][i + 1]->isInside(step * t2 + O, limit_length));
                    }

                    if (condition)
                    {

                        // if the collision are too close or negative.
                        if (Walker::bouncing)
                        {
                            if (t2 >= EPS_VAL && t2 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t2);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                        else
                        {
                            if (t2 >= 0 && t2 <= step_lenght + barrier_tickness)
                            {
                                dist_intersections.push_back(t2);
                                sph_ids.push_back(i);
                                branch_ids.push_back(list_id);
                                cs.push_back(c);
                            }
                        }
                    }
                }
            }
        }
    }
    

    if (dist_intersections.size() > 0)
    {
        auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
        unsigned index_ = std::distance(std::begin(dist_intersections), min_distance_int);

        int sphere_ind = sph_ids[index_];
        int branch_ind = branch_ids[index_];

        double dist_to_collision = dist_intersections[index_];

        colision.type = Collision::hit;
        colision.rn = cs[index_];

        if (walker.initial_location == Walker::intra)
        {
            if (check_inside_(all_cs))
                colision.col_location = Collision::inside;
            else
                colision.col_location = Collision::outside;
        }
        else if (walker.initial_location == Walker::extra)
        {
            if (check_outside_(all_cs))
                colision.col_location = Collision::outside;
        
            else
                colision.col_location = Collision::inside;
        }

        colision.t = fmin(dist_to_collision, step_lenght);
        colision.colision_point = walker.pos_v + colision.t * step;

        Vector3d normal = (colision.colision_point - spheres_list[branch_ind][sphere_ind]->center).normalized();

        Vector3d temp_step = step;
        elasticBounceAgainsPlane_intra(walker.pos_v, normal, colision.t, temp_step);
        colision.bounced_direction = temp_step.normalized();

        return true;
    }
    else
    {
        colision.type = Collision::null;
        return false;
    }
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
    cout << "intersection_sphere_vector" << endl;
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

void Neuron::set_dendrites(std::vector<Dendrite> &dendrites_to_add)
{
    if (dendrites_to_add.size() != 0)
        for(size_t i=0; i < dendrites_to_add.size(); ++i)
            add_dendrite(dendrites_to_add[i]);
}

void Neuron::add_dendrite(Dendrite &dendrite_to_add)
{
    dendrite_to_add.add_projection();
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