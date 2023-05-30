#include "mcsimulation.h"
#include "Eigen/Dense"
#include "simerrno.h"
#include "pgsesequence.h"
#include "gradientwaveform.h"
#include "dynamic_Cylinder.h"

using namespace std;

int MCSimulation::count =0;

MCSimulation::MCSimulation()
{
    plyObstacles_list = nullptr;
    dynamicsEngine = nullptr;
    dataSynth = nullptr;
    sphere_list = nullptr;
    cylinder_list = nullptr;
    dyn_cylinder_list = nullptr;
    dyn_sphere_list = nullptr;
    axon_list= nullptr;
    neuron_list= nullptr;

    id = count;
    count++;
}

/*DEPRECATED*/
MCSimulation::MCSimulation(std::string config_file)
{
    plyObstacles_list = nullptr;
    dynamicsEngine = nullptr;
    dataSynth      = nullptr;
    sphere_list    = nullptr;
    cylinder_list = nullptr;
    dyn_cylinder_list = nullptr;
    dyn_sphere_list = nullptr;
    axon_list= nullptr;
    neuron_list= nullptr;


    params.readSchemeFile(config_file);
    dynamicsEngine = new DynamicsSimulation(params);

    if(params.scheme_file.length() > 2){
        scheme.readSchemeFile(params.scheme_file,params.scale_from_stu);
    }


    if(scheme.type == "PGSE"){
        dataSynth = new PGSESequence(scheme);
        dataSynth->setNumberOfSteps(dynamicsEngine->params.num_steps);

        if(params.subdivision_flag){
            dataSynth->subdivision_flag = true;
            dataSynth->subdivisions = params.subdivisions;
            dataSynth->initializeSubdivisionSignals();
        }
    }

    dataSynth->separate_signal = params.separate_signals;

    dynamicsEngine->id = count;
    id = count;
    count++;
}

MCSimulation::MCSimulation(Parameters& params_)
{
    plyObstacles_list = nullptr;
    dynamicsEngine    = nullptr;
    dataSynth         = nullptr;
    sphere_list       = nullptr;
    cylinder_list     = nullptr;
    dyn_cylinder_list = nullptr;
    dyn_sphere_list   = nullptr;
    axon_list         = nullptr;


    params = params_;
    dynamicsEngine = new DynamicsSimulation(params);

    if(params.scheme_file.length() > 2){
        scheme.readSchemeFile(params.scheme_file,params.scale_from_stu);
    }

    if(scheme.type == "PGSE"){
        dataSynth = new PGSESequence(scheme);
    }
    if(scheme.type == "WAVEFORM"){
        dataSynth = new GradientWaveform(scheme);
    }

    dataSynth->setNumberOfSteps(dynamicsEngine->params.num_steps);
    dataSynth->separate_signal = params.separate_signals;

    if(params.subdivision_flag){
        dataSynth->subdivision_flag = true;
        dataSynth->subdivisions = params.subdivisions;
        dataSynth->initializeSubdivisionSignals();
    }

    if(params.separate_signals)
        dataSynth->initializeIntraExtraSignals();

    dynamicsEngine->id = count;
    id = count;
    count++;
}


void MCSimulation::startSimulation()
{

    iniObstacles();

    if(dataSynth != nullptr){
        dynamicsEngine->startSimulation(dataSynth);
    }
    else{
        dynamicsEngine->startSimulation();
    }

}

double MCSimulation::getExpectedFreeeDecay(unsigned i)
{
    if(dataSynth){
        double b = dataSynth->getbValue(i);
        return exp(-b*params.diffusivity);
    }

    return -1;
}


void MCSimulation::iniObstacles()
{

    addObstacles();

    addPLYObstacles();

    addVoxels();
}

void MCSimulation::addObstacles()
{
    this->dynamicsEngine->cylinders_list = this->cylinder_list;
    this->dynamicsEngine->spheres_list   = this->sphere_list;
    this->dynamicsEngine->dyn_cylinders_list = this->dyn_cylinder_list;
    this->dynamicsEngine->dyn_spheres_list = this->dyn_sphere_list;
    this->dynamicsEngine->axons_list = this->axon_list;
    this->dynamicsEngine->neurons_list = this->neuron_list;
}


////* Auxiliare method to split words in a line using the spaces*//
template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


void MCSimulation::addPLYObstacles()
{

    dynamicsEngine->plyObstacles_list = this->plyObstacles_list;
}

void MCSimulation::addVoxels()
{
    for(unsigned i = 0 ; i < params.voxels_list.size(); i++){
        dynamicsEngine->voxels_list.push_back(Voxel(params.voxels_list[i].first,params.voxels_list[i].second));
    }
}

bool cylinderIsCloseBoundery(Cylinder& cyl, Eigen::Vector3d min_limits,Eigen::Vector3d max_limits){

    //gap to the boundary
    double gap = 1e-6;
    //3 dimensional vector
    for (int i = 0 ; i < 3; i++)
        if( (cyl.P[i] - cyl.radius - gap < min_limits[i]) || (cyl.P[i] + cyl.radius + gap  > max_limits[i]) )
            return true;

    return false;
}

bool dyncylinderIsCloseBoundery(Dynamic_Cylinder& cyl, Eigen::Vector3d min_limits,Eigen::Vector3d max_limits){

    //gap to the boundary
    double gap = 1e-6;
    //3 dimensional vector
    for (int i = 0 ; i < 3; i++)
        if( (cyl.P[i] - cyl.radius - gap < min_limits[i]) || (cyl.P[i] + cyl.radius + gap  > max_limits[i]) )
            return true;

    return false;
}


MCSimulation::~MCSimulation()
{
    if(dynamicsEngine != nullptr)
        delete dynamicsEngine;

    if(dataSynth != nullptr)
        delete dataSynth;
}
