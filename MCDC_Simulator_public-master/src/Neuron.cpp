#include "Neuron.h"
#include "dynamic_sphere.h"
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

    ub = 5; // upper bound for span_radius 
    lb = 2; // lower bound for span_radius 
    uniform_int_distribution<mt19937::result_type> dist_span_radius(lb, ub);

    // Generate int number in [lb, ub], in [mm]
    span_radius = dist_span_radius(rng);
    span_radius /=10;

}

Neuron::~Neuron()
{
    nb_neurons--;
}

Neuron::Neuron(vector<Axon> dendrites_, Sphere soma_) : Neuron()
{
    dendrites = dendrites_;          
    soma = soma_;
}

Neuron::Neuron(Vector3d soma_center, double soma_radius=5e-3) : Neuron()
{
    soma = Sphere(soma_center, soma_radius);
}

Neuron::Neuron(vector<Axon> dendrites_, Vector3d soma_center, double soma_radius=5e-3) : Neuron()
{
    dendrites = dendrites_;          
    soma = Sphere(soma_center, soma_radius);
}

Neuron::Neuron(const Neuron &neuron)
{
    id = nb_neurons++;
    dendrites = neuron.dendrites;           
    soma = neuron.soma;
    nb_dendrites = neuron.nb_dendrites;
    span_radius = neuron.span_radius;
}

double Neuron::minDistance(Walker &w)
{
    vector<double> distances;
    distances.clear();
    distances = Distances_to_Spheres(w);
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

bool Neuron::isPosInsideNeuron(Eigen::Vector3d &position,  double barrier_thickness, bool swell_){
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
    tie(neuron_part, part_id) = isNearNeuron(position, barrier_thickness);
    if(!(neuron_part=="none")){
        if (neuron_part=="dendrite")
        {
            // If the neurite is swelling, check collision with the maximum radius
            if (swell_){rad = dendrites[part_id].max_radius;}
            else{rad = dendrites[part_id].radius;}
            
            // vect size 3 of vect of projection_pt
            coliding_projs = dendrites[part_id].projections.find_collisions_all_axes(position, rad, dendrites[part_id].id);
            
            if (coliding_projs.size() == 3){ 

                // for all coliding objects in x 
                for(unsigned j = 0; j < coliding_projs[0].size() ; j++)
                { 
                    const Projections::projection_pt coliding_proj = coliding_projs[0][j];
                    // if the same coliding objects are also in y and z but are not from same objects
                    if (swell_){
                        colliding_all_axes = (dendrites[part_id].projections_max.isProjInside(coliding_projs[1], coliding_proj) && 
                                              dendrites[part_id].projections_max.isProjInside(coliding_projs[2], coliding_proj));
                    }
                    else{
                        colliding_all_axes = (dendrites[part_id].projections.isProjInside(coliding_projs[1], coliding_proj) && 
                                              dendrites[part_id].projections.isProjInside(coliding_projs[2], coliding_proj));
                    }
                    if (colliding_all_axes){
                        sphere_ = dendrites[part_id].spheres[coliding_proj.sph_id];
                        if (swell_)
                        {
                            rad_ = sphere_.max_radius;
                        }
                        else
                        {
                            rad_ = sphere_.radius;
                        }
                        if (sphere_.distSmallerThan(position, -barrier_thickness + rad_))
                        { 
                            return true;
                            break;
                        }  
                    } 
                }
            }
        } 
        else if (neuron_part == "soma")
        {
            double dis = soma.minDistance(position);
            if( dis <= barrier_thickness ){return true;}
        }  
    }
    return false;
} 


tuple<string, int> Neuron::isNearNeuron(Vector3d &position,  double barrier_thickness)
{
    // Check soma box
    for (unsigned int axis=0 ; axis < 3 ; ++axis)
    {
        if (position[axis] >= soma.center[axis] - soma.radius + barrier_thickness && position[axis] <= soma.center[axis] + soma.radius - barrier_thickness)
        {
            return tuple<string, int>{"soma", soma.id}; 
        }
    }
    
    // Check each dendrite's box
    bool isnear = false;
    for (unsigned int i=0 ; i < nb_dendrites ; ++i)
    {
        for (unsigned int axis=0 ; axis < 3 ; ++axis)
        {
            Vector2d axis_limits = dendrites[i].projections.axon_projections[axis];

            if (position[axis] >= axis_limits[0] + barrier_thickness && position[axis] <= axis_limits[1] - barrier_thickness)
            {
                isnear = isnear && true;
            }
        }
        if (isnear)
        {
            return tuple<string, int>{"dendrite", i};
        }
    }

    return tuple<string, int>{"none", -1};
}


bool Neuron::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{

    string message;
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d next_step = step*step_lenght+O;

    bool isintra;
    bool next_step_intra = isPosInsideNeuron(next_step, barrier_tickness, false);
    bool collision_check = false;


    if (walker.location == Walker::intra){
        isintra = true;
    }
    else if (walker.location == Walker::extra) {
        isintra = false;
    }
    else{
        if(walker.initial_location == Walker::intra){
            isintra = true;
        }
        else{
            isintra = false;
        }
    }

    // if is intra and so is next step -> no collision with border
    if(isintra && next_step_intra){
        colision.type = Collision::null;
        return false;
    }

    // distances to intersections
    std::vector<double> dist_intersections;
    std::vector<double> cs;
    std::vector<int> sph_ids;
    sph_ids.clear();
    dist_intersections.clear();
    int sph_id;

    // find a sphere that is near the walker 
    string neuron_part; // "soma", "dendrite" or "none"
    int part_id; // id of the soma or dendrite. -1 if not in neuron
    int closest_sphere_index; // id of the sphere of the dendrite
    tie(neuron_part, part_id, closest_sphere_index) = closest_sphere_dichotomy(walker, step_lenght, barrier_tickness);

    if (neuron_part == "dendrite")
    {
        int sph_begin;
        int sph_end;
        // Check the spheres that are in the vicinity of the closest sphere
        sph_begin = closest_sphere_index - int((step_lenght + barrier_tickness + dendrites[part_id].radius)*5 / dendrites[part_id].radius)-1;
        if(sph_begin < 0){
            sph_begin = 0;
        } 
        sph_end  = closest_sphere_index + int((step_lenght + barrier_tickness + dendrites[part_id].radius)*5 / dendrites[part_id].radius)+1;
        if (sph_end > dendrites[part_id].spheres.size()){
            sph_end = dendrites[part_id].spheres.size();
        } 

        for (unsigned i = sph_begin ; i < sph_end; ++i){
            // distances to collision
            double t1;
            double t2;
            double c;
            bool intersect = intersection_sphere_vector(t1, t2,  dendrites[part_id].spheres[i], step, step_lenght, O, isintra, c);
            if (intersect){
                if(walker.status == Walker::bouncing){
                    //if the collision are too close or negative.
                    if(  t1 >= EPS_VAL && t1 <= step_lenght + barrier_tickness){
                        dist_intersections.push_back(t1);
                        sph_id = i;
                        sph_ids.push_back(sph_id);
                        cs.push_back(c);
                    }
                    if(  t2 >= EPS_VAL && t2 <= step_lenght + barrier_tickness){
                        dist_intersections.push_back(t2);
                        sph_id = i;
                        sph_ids.push_back(sph_id);
                        cs.push_back(c);
                    }
                }
                else{
                    if( t1 >= 0 && t1 <= step_lenght + barrier_tickness){
                        dist_intersections.push_back(t1);
                        sph_id = i;
                        sph_ids.push_back(sph_id);
                        cs.push_back(c);
                    }
                    if( t2 >= 0 && t2 <= step_lenght + barrier_tickness){
                        dist_intersections.push_back(t2);
                        sph_id = i;
                        sph_ids.push_back(sph_id);
                        cs.push_back(c);
                        
                    }
                }
            }
        }

        if(dist_intersections.size() > 0){

            unsigned index_ ;

            if (!isintra){
                auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
                index_ = std::distance(std::begin(dist_intersections), min_distance_int);
            }
            else{
                auto max_distance_int = std::max_element(std::begin(dist_intersections), std::end(dist_intersections));
                index_ = std::distance(std::begin(dist_intersections), max_distance_int);
            }

            int sphere_ind = sph_ids[index_];

            double dist_to_collision = dist_intersections[index_];

            Dynamic_Sphere colliding_sphere = dendrites[part_id].spheres[sphere_ind];

            colision.type = Collision::hit;
            colision.rn = cs[index_];

            if(colision.rn <-1e-10)
            {
                colision.col_location = Collision::inside;
                walker.in_obj_index = -1;
            }
            else if(colision.rn >1e-10)
            { 
                colision.col_location = Collision::outside;
            }
            else
            {
                colision.col_location = Collision::unknown;
            }
        
            colision.t = fmin(dist_to_collision,step_lenght);
            colision.colision_point = walker.pos_v + colision.t*step;
                
            //Normal point
            Vector3d normal = (colision.colision_point- colliding_sphere.center).normalized();
            Vector3d temp_step = step;
            elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);
            if ((isintra && colision.rn >= -1e-10)){
                colision.col_location = Collision::inside;
                colision.bounced_direction = (-step).normalized();
            }

            else{  
                colision.bounced_direction = temp_step.normalized();
            } 
            return true;
            
        }
        else{
            colision.type = Collision::null;
            return false;   
        }
    }
    else if (neuron_part == "soma")
    {
        
    }
}

bool Neuron::intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere &s, Vector3d &step, double &step_length, Vector3d &pos, bool isintra, double &c){
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    Vector3d m = pos - s.center;
    double rad = s.radius;

    // collision distance
    double d_ = m.norm() - rad;

    //If the minimum distance from the walker to the cylinder is more than
    // the actual step size, we can discard this collision.
    if(d_> EPS_VAL){
        if(d_ > step_length+barrier_tickness){
            return false;
        }
    }

    double a = 1;
    double b = (m.dot(step));
    c = m.dot(m) - rad*rad;

    
    double discr = b*b - a*c;

    if (discr < 0.0 ){
        return false;
    }

    t1 = (-b + sqrt(discr))/(a);
    t2 = (-b - sqrt(discr))/(a);

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