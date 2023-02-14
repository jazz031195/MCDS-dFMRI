#include "neurondistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;

NeuronDistribution::NeuronDistribution(unsigned num_obstacles_, double icvf_, Eigen::Vector3d & min_limits_, Eigen::Vector3d & max_limits_)
{
    num_obstacles = num_obstacles_;
    icvf = icvf_;
    min_limits = min_limits_;
    max_limits = max_limits_;
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
    double best_icvf = 0;
    Eigen::Vector3d best_max_limits;
    min_limits = {0.,0.,0.};

    bool achieved = false;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);

    int tried = 0;

    uint adjustments = 0;
     // We increease 1% the total area. (Is prefered to fit all the spheres than achieve a perfect ICVF.)
    double adj_increase = icvf*0.01;
    while(!achieved){

        double target_icvf = this->icvf+adjustments*adj_increase;

        for(uint t = 0 ;  t < repetition; t++){
            vector<Neuron> neurons_to_add;

            neurons.clear();
            for(unsigned i = 0 ; i < num_obstacles; i++){
                unsigned stuck = 0;
                while(++stuck <= 1000){
                    achieved = true;
                    double t = udist(gen);
                    double x = (t*max_limits[0] + (1-t)*min_limits[0]);
                    t = udist(gen);
                    double y = (t*max_limits[1] + (1-t)*min_limits[1]);
                    t = udist(gen);
                    double z = (t*max_limits[2] + (1-t)*min_limits[2]);

                    Eigen::Vector3d soma_center = {x, y, z};
                    double soma_radius = 5e-3; //mm
                    Neuron neuron(soma_center, soma_radius);
                    neurons_to_add.push_back(neuron);

                    double min_distance;

                    bool collision = false;
                    //checkForCollition(neuron, min_limits, max_limits, neurons_to_add, min_distance);

                    if(!collision){
                        for (unsigned j = 0; j < neurons_to_add.size(); j++){
                            neurons.push_back(neurons_to_add[j]);
                        }
                        break;
                    }
                    
                }

                // int dummy;
                // double icvf_current = computeICVF(spheres,min_limits, max_limits,dummy);
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
    // double icvf_current = computeICVF(spheres,min_limits, max_limits,perc_);
    // string  message = "Percentage of spheres  selected: "+ to_string(double(perc_)/radiis.size()*100.0)
    //        + "%,\nICVF achieved: " + to_string(icvf_current*100) + "  ("+ to_string( int((icvf_current/icvf*100))) + "% of the desired icvf)\n";
    // SimErrno::info(message,cout);
}

void NeuronDistribution::printSubstrate(ostream &out)
{
    out << 1e-3 << endl;
    out << icvf << endl;

    for (unsigned i = 0; i < neurons.size(); i++)
    {
        SimErrno::info(std::to_string(neurons[i].soma.center[0]),cout);

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

std::vector<NeuronDistribution::projection_pt> NeuronDistribution::find_collisions(projection_pt proj_on_axis_min, projection_pt proj_on_axis_max,std::vector<projection_pt> projections_on_axis, ostream& out)
{   
   
}

bool NeuronDistribution::search_for_sphere(std::vector<NeuronDistribution::projection_pt> spheres_, NeuronDistribution::projection_pt s){
}  

bool NeuronDistribution::isSphereColliding(Dynamic_Sphere sph, double distance_to_be_inside, int axon_id, int sph_id, ostream& out)
{
} 
 
bool NeuronDistribution::isColliding(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out){
   
}


double NeuronDistribution::computeICVF(std::vector<Axon> &axons, Vector3d &min_limits, Vector3d &max_limits, int &num_no_repeat)
{

}



// std::tuple<double, double>  phi_gamma_to_target(Eigen::Vector3d prev_pos, Eigen::Vector3d new_pos, Eigen::Vector3d end,  ostream& out) 
// {

// }

bool NeuronDistribution::check_borders(Eigen::Vector3d pos, double distance_to_border) 
{
    
}

std::vector<Dynamic_Sphere> NeuronDistribution::GrowAxon(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out)
{

}  