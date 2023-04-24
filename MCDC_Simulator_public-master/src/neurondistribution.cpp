#include "neurondistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;

NeuronDistribution::NeuronDistribution(int const& num_obstacles_, double const& icvf_, Vector3d const& min_limits_vx_, Vector3d const& max_limits_vx_, double const& step_length_):
min_limits_vx(min_limits_vx_), max_limits_vx(max_limits_vx_), num_obstacles(num_obstacles_), icvf(icvf_), step_length(step_length_)
{
    neurons.clear();
    projections_x.clear();
    projections_y.clear();
    projections_z.clear();
    
    string message = "neurons : " + std::to_string(this->num_obstacles) + " \n";
    SimErrno::info(message, std::cout);
}

// void NeuronDistribution::computeMinimalSize(std::vector<double> const& radiis, double &icvf_, Vector3d &l) const
// {

//     /*A little heuristic for complicated ICVF: > 0.7*/
//     if (icvf_ >= 0.7 && icvf_ < 0.99)
//         icvf_ += 0.01;

//     double area = 0;

//     for (uint i = 0; i < radiis.size(); i++)
//         area += radiis[i] * radiis[i] * M_PI;

//     double l_ = sqrt(area / icvf_);

//     l = {l_, l_, l_};
// }


void NeuronDistribution::createSubstrate()
{

    uint repetition = 1;
    uint max_adjustments = 5;
    // double best_icvf = 0;
    // Eigen::Vector3d best_max_limits;

    bool achieved = false;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);

    uint adjustments = 0;
     // We increease 1% the total area. (Is prefered to fit all the spheres than achieve a perfect ICVF.)
    // double adj_increase = icvf*0.01;
    while(!achieved){

        // double target_icvf = this->icvf+adjustments*adj_increase;
        double soma_radius = 10e-3; //mm
        // Let enough distance for the radius and for a step_length so that 
        // mirroring border conditions are ok
        double min_distance_from_border = barrier_tickness + soma_radius;
        for(uint t = 0 ;  t < repetition; t++){
            neurons.clear();
            vector<Eigen::Vector3d> soma_centers;
            soma_centers.clear();
            
            for(int i = 0 ; i < num_obstacles; i++){
                unsigned stuck = 0;
                while(++stuck <= 10000){
                    double t = udist(gen);
                    double x = (t*max_limits_vx[0] + (1-t)*min_limits_vx[0]);
                    t        = udist(gen);
                    double y = (t*max_limits_vx[1] + (1-t)*min_limits_vx[1]);
                    t        = udist(gen);
                    double z = (t*max_limits_vx[2] + (1-t)*min_limits_vx[2]);

                    Eigen::Vector3d soma_center = {x, y, z};
                   
                    // If too close to the border, discard
                    if ((x - min_limits_vx[0] <= min_distance_from_border) ||
                        (max_limits_vx[0] - x <= min_distance_from_border) ||
                        (y - min_limits_vx[1] <= min_distance_from_border) ||
                        (max_limits_vx[1] - y <= min_distance_from_border) ||
                        (z - min_limits_vx[2] <= min_distance_from_border) ||
                        (max_limits_vx[2] - z <= min_distance_from_border))
                        {
                            continue;
                        } 
                    if (!isSphereColliding(soma_center, soma_radius, soma_centers))
                    {
                        soma_centers.push_back(soma_center);
                        break;
                    }  
                }
            } // end for neurons
            for(size_t i=0; i < soma_centers.size(); ++i)
            {
                Neuron neuron(soma_centers[i], soma_radius);
                growDendrites(neuron); 
                neurons.push_back(neuron);
                cout << " End of neuron " << i << endl;
            }

            double icvf, somaFraction, dendritesFraction;
            tie(icvf, somaFraction, dendritesFraction) = computeICVF(0);
            cout << icvf << somaFraction << dendritesFraction << endl;
            achieved = true;

            // if(this->icvf - best_icvf  < 0.0005){
            //     achieved = true;
            //     break;
            // }
        }

        adjustments++;
        if(adjustments > max_adjustments){
            break;
        }
    }

    // spheres = best_spheres;
    // max_limits = best_max_limits;

    // //TODO cambiar a INFO
    // int perc_;
    // double icvf_current = computeICVF(spheres,min_limits_vx, max_limits_vx,perc_);
    // string  message = "Percentage of spheres  selected: "+ to_string(double(perc_)/radiis.size()*100.0)
    //        + "%,\nICVF achieved: " + to_string(icvf_current*100) + "  ("+ to_string( int((icvf_current/icvf*100))) + "% of the desired icvf)\n";
    // SimErrno::info(message,cout);
}

void NeuronDistribution::growDendrites(Neuron& neuron)
{
    // Store all the starting points of dendrites, on the soma of the neuron
    std::vector<Eigen::Vector3d> start_dendrites;
    int max_tries = 10000;

    for(uint8_t i = 0; i < neuron.nb_dendrites; ++i)
    {   
        int tries = 0;
        int nb_branching = generateNbBranching();
        // Radius of each dendrite sphere [mm]
        double sphere_radius = 0.5e-3;
        // Don't initiate dendrite too close from the borders
        double min_distance_from_border = barrier_tickness + sphere_radius;
        
        while(tries < max_tries)
        {
            Vector3d dendrite_start = generatePointOnSphere(neuron.soma.center, neuron.soma.radius);

            while(!isInVoxel(dendrite_start, min_distance_from_border) || (isSphereColliding(dendrite_start, sphere_radius)))
            {
               dendrite_start = generatePointOnSphere(neuron.soma.center, neuron.soma.radius);
            }

            // If the vector is not already contained in start_dendrites, add it. 
            // Otherwise, decrement i and do one more round
            if((i != 0) && 
               std::count(start_dendrites.begin(), start_dendrites.end(), dendrite_start)){ i--; tries++; }
            else
            {
                start_dendrites.push_back(dendrite_start);
                Eigen::Vector3d dendrite_direction = dendrite_start - neuron.soma.center;
                dendrite_direction.normalize();
                Dendrite dendrite;
                // Tuple xyz center and center id
                vector<tuple<Vector3d, int>> parent_centers {{dendrite_start, 0}};
                
                // Create the subbranches
                for(int b=0; b < nb_branching; ++b)
                {
                    // Length of a segment before branching
                    int l_segment = generateLengthSegment();
                    // Number of spheres per segment
                    int nb_spheres = l_segment / (sphere_radius/4); //Let's assume that dendrites have a radius of 0.5microns so far

                    if(b == 0)
                    {
                        vector<int> proximal_branch {0};
                        vector<int> distal_branch {1, 2};
                        growSubbranch(dendrite, parent_centers[0], dendrite_direction, nb_spheres, sphere_radius, proximal_branch, distal_branch, 
                                      min_distance_from_border);
                    }

                    for(size_t p=0; parent_centers.size(); p++)
                    {
                        Eigen::Vector3d begin;
                        Axon subbranch(sphere_radius, begin, begin, 0, false, false , 1);
                        std::vector<Dynamic_Sphere> spheres_to_add;
                        spheres_to_add.clear();

                        Eigen::Vector3d center = {0, 0, 0};
                        bool discard_dendrite  = false;
                        
                        Vector3d origin;
                        int parent_id;
                        tie(origin, parent_id) = parent_centers[p];

                        for(int j=0; j < nb_spheres; ++j)
                        {
                            center = j*dendrite_direction*sphere_radius/4 + origin;

                            if(isInVoxel(center, min_distance_from_border))
                            {
                                if (!isSphereColliding(center, sphere_radius))
                                {
                                    Dynamic_Sphere sphere_to_add(center, sphere_radius, 0, false, j, 1, false);
                                    spheres_to_add.push_back(sphere_to_add);
                                }
                            }
                            // else
                            // {
                            //     discard_dendrite = true;
                            //     createTwinSphere(center, sphere_radius, discard_dendrite, j);
                            //     if(!discard_dendrite)
                            //     {
                            //         Dynamic_Sphere sphere_to_add(center, sphere_radius, 0, false, j, 1, false);
                            //         spheres_to_add.push_back(sphere_to_add);
                            //     }
                            //     else
                            //         break;
                            // } 
                        }
                        if (spheres_to_add.size() == 0)
                            break;
                        if (!discard_dendrite)
                        {
                            subbranch.set_spheres(spheres_to_add, i);
                            dendrite.add_subbranch(subbranch);
                            neuron.add_dendrite(dendrite);
                            break;
                        }  
                    }   
                }  
            }
        }
    } 
}

void NeuronDistribution::growSubbranch(Dendrite& dendrite, tuple<Vector3d, int> const& parent, Vector3d const& dendrite_direction, 
                                      int const& nb_spheres, double const& sphere_radius, vector<int> const& proximal_end, 
                                      vector<int> const& distal_end, double const& min_distance_from_border)
{
    Eigen::Vector3d begin;
    Axon subbranch(sphere_radius, begin, begin, 0, false, false , 1);
    std::vector<Dynamic_Sphere> spheres_to_add;
    spheres_to_add.clear();

    Eigen::Vector3d center = {0, 0, 0};
    bool discard_dendrite  = false;
    
    Vector3d origin_branch;
    int parent_id;
    tie(origin_branch, parent_id) = parent;

    for(int j=0; j < nb_spheres; ++j)
    {
        center = j*dendrite_direction*sphere_radius/4 + origin_branch;

        if(isInVoxel(center, min_distance_from_border))
        {
            if (!isSphereColliding(center, sphere_radius))
            {
                Dynamic_Sphere sphere_to_add(center, sphere_radius, 0, false, j, 1, false);
                spheres_to_add.push_back(sphere_to_add);
            }
        }
    }

    if (!discard_dendrite)
    {
        int subbranch_id = dendrite.get_nb_subbranches() + 1;
        subbranch.set_spheres(spheres_to_add, subbranch_id);
        dendrite.add_subbranch(subbranch);
    }  
}

int NeuronDistribution::generateNbBranching(int const& lower_bound, int const& upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<int> nbBranching(lower_bound, upper_bound);

    return nbBranching(rng);
}

int NeuronDistribution::generateLengthSegment(double const& lower_bound, double const& upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<double> segmentLength(lower_bound, upper_bound);

    return segmentLength(rng);
}

Vector3d NeuronDistribution::generatePointOnSphere(Vector3d const& center, double const& radius) const
{
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::normal_distribution<double> distribution(0.0, 1.0);

    // Find a point at the surface of the soma
    double x = distribution(generator); 
    double y = distribution(generator);
    double z = distribution(generator);

    // Avoid division by 0 while normalizing
    while(x==0 && y==0 && z==0)
    {
        x = distribution(generator);
        y = distribution(generator);
        z = distribution(generator);
    }
    double normalization_factor = sqrt(x*x + y*y + z*z);
    x = x/normalization_factor*radius + center[0];
    y = y/normalization_factor*radius + center[1];
    z = z/normalization_factor*radius + center[2];
    
    return {x, y, z};   
}

void NeuronDistribution::createTwinSphere(Vector3d &center, double const& sphere_radius, bool &discard_dendrite, size_t const& j)
{
    Vector3d new_center = center;
    for (size_t axis=0; axis < 3; ++axis)
    {
        double distance_from_min_lim = center[axis] - min_limits_vx[axis];
        if (distance_from_min_lim < sphere_radius)
        {
            new_center[axis] = max_limits_vx[axis] + distance_from_min_lim;
        }
        double distance_from_max_lim = max_limits_vx[axis] - center[axis];
        if (distance_from_max_lim < sphere_radius)
        {
            new_center[axis] = min_limits_vx[axis] - distance_from_max_lim;
        }
    }
    if(!isSphereColliding(new_center, sphere_radius))
    {
        discard_dendrite = false;
        center = new_center;
    }
}

bool NeuronDistribution::isSphereColliding(Dynamic_Sphere const& sph) 
{
    Vector3d position = sph.center;
    double distance_to_be_inside = sph.max_radius + 2 * barrier_tickness;
    int dummy, dummy2;
    for (unsigned i = 0; i < neurons.size() ; i++){
        bool isinside = neurons[i].isPosInsideNeuron(position, distance_to_be_inside, false, dummy, dummy2);
        if (isinside)
            return true;
    }
    return false;
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius, vector<Vector3d> const& soma_centers) 
{
    double distance_to_be_inside = 2 * barrier_tickness;
    for (unsigned i = 0; i < soma_centers.size() ; i++){
        bool isinside = isSphereCollidingSphere(sphere_center, soma_centers[i], sphere_radius, sphere_radius, distance_to_be_inside);
        if (isinside)
            return true;
    }
    return false;
}

bool NeuronDistribution::isSphereCollidingSphere(Vector3d const& pos1, Vector3d const& pos2, double const& radius1, double const& radius2, double const& minDistance) const 
{
    Vector3d m = pos1 - pos2;
    double distance_to_sphere = m.norm() - radius1 - radius2;

    return distance_to_sphere < minDistance;
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius) 
{
    double distance_to_be_inside = sphere_radius + 2 * barrier_tickness;
    int dummy, dummy2;
    for (unsigned i = 0; i < neurons.size() ; i++){
        bool isinside = neurons[i].isPosInsideNeuron(sphere_center, distance_to_be_inside, false, dummy, dummy2);
        if (isinside)
            return true;
    }
    return false;
}

void NeuronDistribution::printSubstrate(ostream &out) const
{
    out << 1 << endl; //scale
    out << 0 << endl; //volume_inc_perc
    out << 0 << endl; //dyn_perc
    out << icvf << endl;
    out << min_limits_vx[0] << endl; //min_limits [mm]
    out << max_limits_vx[0] << endl; //max_limits [mm]

    for (unsigned i = 0; i < neurons.size(); i++)
    {
        // Print for soma : x y z r bool_active
        // bool_active = 1 if the sphere can be activated (swollen)
        out << neurons[i].soma.center[0] << " " 
        << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius    <<  " "
        << 0 << endl; //bool_active = false for now
        out << "Soma " + to_string(i) << endl;

        for (size_t j = 0; j < neurons[i].dendrites.size(); j++)
        {
            for (size_t k = 0; k < neurons[i].dendrites[j].subbranches.size(); k++)
            {
                for (size_t l = 0; l < neurons[i].dendrites[j].subbranches[k].spheres.size(); l++)
                {
                    // Print for each dendrite, each sphere
                    out << neurons[i].dendrites[j].subbranches[k].spheres[l].center[0] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[1] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[2] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].radius << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].swell << endl; 
                }
                out << "Segment " + to_string(k) << endl;
            }
            out << "Dendrite " + to_string(j) << endl;
        }
        out << "Neuron " + to_string(i) << endl;
    }
}

bool NeuronDistribution::isInVoxel(Eigen::Vector3d const& pos, double const& distance_to_border) const
{
    Eigen::Vector3d new_min_limits_vx = {min_limits_vx[0] + distance_to_border, min_limits_vx[1] + distance_to_border, min_limits_vx[2] + distance_to_border};
    Eigen::Vector3d new_max_limits_vx = {max_limits_vx[0] - distance_to_border, max_limits_vx[1] - distance_to_border, max_limits_vx[2] - distance_to_border};

    if ((pos[0] - new_min_limits_vx[0]) < 0 || 
        (pos[1] - new_min_limits_vx[1]) < 0 || 
        (pos[2] - new_min_limits_vx[2]) < 0)
        return false;
    else if ((pos[0] - new_max_limits_vx[0]) > 0 || 
             (pos[1] - new_max_limits_vx[1]) > 0 || 
             (pos[2] - new_max_limits_vx[2]) > 0) 
        return false;

    return true;   
}

tuple<double, double, double> NeuronDistribution::computeICVF(double const& min_distance_from_border) const
{

    if (neurons.size() == 0)
        return make_tuple(0, 0, 0);

    double VolumeVoxel = (max_limits_vx[0] - min_limits_vx[0] - min_distance_from_border) * (max_limits_vx[1] - min_limits_vx[1] - min_distance_from_border) * (max_limits_vx[2] - min_limits_vx[2] - min_distance_from_border);
    double VolumeSoma = 0;
    double VolumeDendrites = 0;
   
    for (size_t i = 0; i < neurons.size(); i++)
    {
        // Calculate the volume of the soma
        VolumeSoma += 4/3*M_PI*pow(neurons[i].soma.radius, 3);

        // Calculate the cylindrical volume of each dendrite
        for (uint8_t j = 0; j < neurons[i].nb_dendrites; j++)
        {
            VolumeDendrites += neurons[i].dendrites[j].volumeDendrite();
        }      
    }
    
    double somaFraction      = VolumeSoma / VolumeVoxel;
    double dendritesFraction = VolumeDendrites/ VolumeVoxel;
    double ICVF              = somaFraction + dendritesFraction;
    return make_tuple(ICVF, somaFraction, dendritesFraction);
}
