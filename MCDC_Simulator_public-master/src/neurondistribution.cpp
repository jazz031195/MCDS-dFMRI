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
        double soma_radius = 5e-3; //mm
        // Let enough distance for the radius and for a step_length so that 
        // mirroring border conditions are ok
        double min_distance_from_border = barrier_tickness + soma_radius + step_length;
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
                    if (!isSphereColliding(soma_center, soma_radius))
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
                
            icvf = computeICVF(min_distance_from_border);
            cout << icvf << endl;
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
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::normal_distribution<double> distribution(0.0, 1.0);
    int max_tries = 10000;

    for(uint8_t i = 0; i < neuron.nb_dendrites; ++i)
    {   
        int tries = 0;

        while(tries < max_tries)
        {
            // Radius of each dendrite sphere [mm]
            double sphere_radius = 0.5e-3;
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
            x = x/normalization_factor*neuron.soma.radius + neuron.soma.center[0];
            y = y/normalization_factor*neuron.soma.radius + neuron.soma.center[1];
            z = z/normalization_factor*neuron.soma.radius + neuron.soma.center[2];
            Eigen::Vector3d dendrite_start(x, y, z);
            
            // Keep the neuron one step length away from the border so that water mirroring at
            // the boundary doesn't need to be checked if intra
            double min_distance_from_border = barrier_tickness + sphere_radius + step_length;

            while(!isInVoxel(dendrite_start, min_distance_from_border) || (isSphereColliding(dendrite_start, sphere_radius)))
            {
                x = distribution(generator);
                y = distribution(generator);
                z = distribution(generator);
                dendrite_start = {x, y, z};
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
                int nb_spheres = neuron.span_radius / (sphere_radius/4); //Let's assume that dendrites have a radius of 0.5microns so far
                
                Eigen::Vector3d begin;
                Axon dendrite(sphere_radius, begin, begin, 0, false, false , 1);
                std::vector<Dynamic_Sphere> spheres_to_add;
                spheres_to_add.clear();

                for(int j=0; j < nb_spheres; ++j)
                {
                    Eigen::Vector3d center = j*dendrite_direction*sphere_radius/4 + dendrite_start;
                    if (isInVoxel(center, min_distance_from_border) && !isSphereColliding(center, sphere_radius))
                    {
                        Dynamic_Sphere sphere_to_add(center, sphere_radius, 0, false, j, 1, false);
                        spheres_to_add.push_back(sphere_to_add);
                    }
                    else
                        break; 
                }
                if (spheres_to_add.size() == 0)
                    break;
                dendrite.set_spheres(spheres_to_add, i);
                neuron.add_dendrite(dendrite);
                break;
            }
        }
    } 
}

bool NeuronDistribution::isSphereColliding(Dynamic_Sphere const& sph) 
{
    Vector3d position = sph.center;
    double distance_to_be_inside = sph.max_radius + 10 * barrier_tickness;
    int dummy, dummy2;
    for (unsigned i = 0; i < neurons.size() ; i++){
        bool isinside = neurons[i].isPosInsideNeuron(position, distance_to_be_inside, false, dummy, dummy2);
        if (isinside)
            return true;
    }
    return false;
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius) 
{
    double distance_to_be_inside = sphere_radius + 10 * barrier_tickness;
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
            for (size_t k = 0; k < neurons[i].dendrites[j].spheres.size(); k++)
            {
                // Print for each dendrite, each sphere
                out << neurons[i].dendrites[j].spheres[k].center[0] << " "
                << neurons[i].dendrites[j].spheres[k].center[1] << " "
                << neurons[i].dendrites[j].spheres[k].center[2] << " "
                << neurons[i].dendrites[j].spheres[k].radius << " "
                << neurons[i].dendrites[j].spheres[k].swell << endl; 
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

double NeuronDistribution::computeICVF(double const& min_distance_from_border) const
{

    if (neurons.size() == 0)
        return 0;

    double VolumeVoxel = (max_limits_vx[0] - min_limits_vx[0] - min_distance_from_border) * (max_limits_vx[1] - min_limits_vx[1] - min_distance_from_border) * (max_limits_vx[2] - min_limits_vx[2] - min_distance_from_border);
    double VolumeNeurons = 0;
   
    for (size_t i = 0; i < neurons.size(); i++)
    {
        // Calculate the volume of the soma
        VolumeNeurons += 4/3*M_PI*pow(neurons[i].soma.radius, 3);

        // Calculate the cylindrical volume of each dendrite
        for (uint8_t j = 0; j < neurons[i].nb_dendrites; j++)
        {
            VolumeNeurons += neurons[i].dendrites[j].volumeAxon();
        }      
    }
    return VolumeNeurons / VolumeVoxel;
}
