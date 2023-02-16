#include "neurondistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;

NeuronDistribution::NeuronDistribution(unsigned num_obstacles_, double icvf_, Eigen::Vector3d & min_limits_vx_, Eigen::Vector3d & max_limits_vx_)
{
    num_obstacles = num_obstacles_;
    icvf = icvf_;
    min_limits_vx = min_limits_vx_;
    max_limits_vx = max_limits_vx_;
    neurons.clear();
    projections_x.clear();
    projections_y.clear();
    projections_z.clear();
    //for checking collision
    
    string message = "neurons : " + std::to_string(this->num_obstacles) + " \n";
    SimErrno::info(message, std::cout);

}

void NeuronDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d &l)
{

    /*A little heuristic for complicated ICVF: > 0.7*/
    if (icvf_ >= 0.7 && icvf_ < 0.99)
    {
        icvf_ += 0.01;
    }

    double area = 0;

    for (uint i = 0; i < radiis.size(); i++)
    {
        area += radiis[i] * radiis[i] * M_PI;
    }

    double l_ = sqrt(area / icvf_);

    l = {l_, l_, l_};
}


void NeuronDistribution::add_projection(Axon ax, int ax_index, double distance_to_be_inside, ostream& out)
{
  
} 

void NeuronDistribution::createSubstrate()
{

    uint repetition = 1;
    uint max_adjustments = 5;
    // double best_icvf = 0;
    // Eigen::Vector3d best_max_limits;
    min_limits_vx = {0.,0.,0.};

    bool achieved = false;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);

    uint adjustments = 0;
     // We increease 1% the total area. (Is prefered to fit all the spheres than achieve a perfect ICVF.)
    double adj_increase = icvf*0.01;
    while(!achieved){

        // double target_icvf = this->icvf+adjustments*adj_increase;

        for(uint t = 0 ;  t < repetition; t++){
            neurons.clear();
            for(unsigned i = 0 ; i < num_obstacles; i++){
                unsigned stuck = 0;
                while(++stuck <= 1000){
                    double t = udist(gen);
                    double x = (t*max_limits_vx[0] + (1-t)*min_limits_vx[0]);
                    t        = udist(gen);
                    double y = (t*max_limits_vx[1] + (1-t)*min_limits_vx[1]);
                    t        = udist(gen);
                    double z = (t*max_limits_vx[2] + (1-t)*min_limits_vx[2]);

                    Eigen::Vector3d soma_center = {x, y, z};
                    double soma_radius = 5e-3; //mm
                    Neuron neuron(soma_center, soma_radius);
                    growDendrites(neuron, i);

                    bool collision = false;
                    //checkForCollition(neuron, min_limits_vx, max_limits_vx, neurons_to_add, min_distance);

                    if(!collision){
                        neurons.push_back(neuron);
                        achieved = true;
                        break;
                    }
                    
                }

                // int dummy;
                // double icvf_current = computeICVF(spheres,min_limits_vx, max_limits_vx,dummy);
                // if(icvf_current > best_icvf ){
                //     best_icvf = icvf_current;
                //     best_spheres.clear();
                //     best_spheres = spheres;
                //     best_max_limits = max_limits;
                // }
            } // end for spheres

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

void NeuronDistribution::growDendrites(Neuron& neuron, int neuron_id)
{
    // Store all the starting points of dendrites, on the soma of the neuron
    std::vector<Eigen::Vector3d> start_dendrites;
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::normal_distribution<double> distribution(0.0, 1.0);
    int max_tries = 10000;

    for(int i = 0; i < neuron.nb_dendrites; ++i)
    {   
        int tries = 0;

        while(tries < max_tries)
        {
            double x = distribution(generator); 
            double y = distribution(generator);
            double z = distribution(generator);

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
            // Radius of each dendrite sphere
            double sphere_radius = 0.5e-3;
            // If the vector is not already contained in start_dendrites, add it. 
            // Otherwise, decrement i and do one more round
            if((i != 0) && std::count(start_dendrites.begin(), start_dendrites.end(), dendrite_start) && (!isInVoxel(dendrite_start, sphere_radius + barrier_tickness))){ i--; tries++; }
            else
            {
                start_dendrites.push_back(dendrite_start);
                Eigen::Vector3d dendrite_direction = dendrite_start - neuron.soma.center;
                dendrite_direction.normalize();
                int nb_spheres = neuron.span_radius / (4*sphere_radius); //Let's assume that dendrites have a radius of 0.5microns so far
                
                Eigen::Vector3d begin;
                Axon dendrite(sphere_radius, begin, begin, 0, false, false , 1);
                std::vector<Dynamic_Sphere> spheres_to_add;

                for(int j=0; j < nb_spheres; ++j)
                {
                    Eigen::Vector3d center = j*dendrite_direction*sphere_radius/4 + dendrite_start;
                    if (isInVoxel(center, sphere_radius + barrier_tickness))
                    {
                        Dynamic_Sphere sphere_to_add(center, sphere_radius, 0, false, j, 1, 0);
                        spheres_to_add.push_back(sphere_to_add);
                    }
                    else{ break; }
                }
                dendrite.set_spheres(spheres_to_add, i);
                dendrite.add_projection(i);
                neuron.dendrites.push_back(dendrite);
                break;
            }
        }
    } 
}

void NeuronDistribution::printSubstrate(ostream &out)
{
    out << 1e-3 << endl;
    out << icvf << endl;

    for (unsigned i = 0; i < neurons.size(); i++)
    {
        out << neurons[i].soma.center[0] << " " << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius << endl;

        for (unsigned j = 0; j < neurons[i].dendrites.size(); j++)
        {
            for (int k = 0; k < neurons[i].dendrites[j].spheres.size(); k++)
            {
                out << neurons[i].dendrites[j].spheres[k].center[0] << " "
                << neurons[i].dendrites[j].spheres[k].center[1] << " "
                << neurons[i].dendrites[j].spheres[k].center[2] << " "
                << neurons[i].dendrites[j].spheres[k].radius << " "
                << endl;
            }
        }
    }
}

bool NeuronDistribution::isInVoxel(Eigen::Vector3d pos, double distance_to_border) 
{
    Eigen::Vector3d new_min_limits_vx = {min_limits_vx[0] + distance_to_border, min_limits_vx[1] + distance_to_border, min_limits_vx[2] + distance_to_border};
    Eigen::Vector3d new_max_limits_vx = {max_limits_vx[0] - distance_to_border, max_limits_vx[1] - distance_to_border, max_limits_vx[2] - distance_to_border};

    if ((pos[0]-new_min_limits_vx[0])<0 || (pos[1]-new_min_limits_vx[1])<0 ){
        return false;
    }
    else if ((pos[0]-new_max_limits_vx[0])>0 || (pos[1]-new_max_limits_vx[1])>0 ) {
        return false;
    }
    return true;   
}