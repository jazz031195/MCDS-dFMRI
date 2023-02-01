#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include "simerrno.h"

using namespace Eigen;
using namespace std;

int Axon::count = 0;
Axon::Axon()
{
    id = count++;
}

Axon::~Axon()
{
    count--;
}

Axon::Axon(const Axon &ax)
{

    swell = ax.swell;
    id = count++;
    spheres = ax.spheres;
    radius = ax.radius;
    begin = ax.begin;
    end = ax.end;
    max_radius = ax.max_radius;
    min_radius = ax.min_radius;
    active_state = ax.active_state;
    projections = ax.projections;
    projections_max = ax.projections_max;
}

void Axon::add_projection(int axon_id){
    Vector2d smallest_pos;
    Vector2d largest_pos;

    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {
        double smallest_pos_ = 100;
        double largest_pos_ = 0;

        for (int i=0; i<spheres.size(); ++i){ 

            double position1, position2, position1_, position2_;
            int sph_id = i;

            // projections_max
            // center + radius
            position1 = spheres[i].center[axis] + max_radius;
            Projections::projection_pt p1 {position1, axon_id, sph_id};
            if (position1 < smallest_pos_){
                smallest_pos_ = position1;
            }
            if (position1 > largest_pos_){
                largest_pos_= position1;
            }
            // center - radius
            position2 = spheres[i].center[axis] - max_radius;
            Projections::projection_pt p2 {position2, axon_id, sph_id};

            if (position2 < smallest_pos_){
                smallest_pos_ = position2;
            }
            if (position2 > largest_pos_){
                largest_pos_= position2;
            }
            projections_max.append_right_place(p1,  p2, axis);

            // projections
            // center + radius
            position1_ = spheres[i].center[axis] + radius;
            Projections::projection_pt p1_ {position1, axon_id, sph_id};
            // center - radius
            position2_ = spheres[i].center[axis] - max_radius;
            Projections::projection_pt p2_ {position2, axon_id, sph_id};

            projections.append_right_place(p1_,  p2_, axis);
        }

        if (axis < 2){
            smallest_pos[axis] = smallest_pos_;
            largest_pos[axis] = largest_pos_;
        }
    }

    Vector2d x_limits = {smallest_pos[0], largest_pos[0]};
    Vector2d y_limits = {smallest_pos[1], largest_pos[1]};
    projections_max.axon_projections.push_back(x_limits);
    projections_max.axon_projections.push_back(y_limits);
    projections.axon_projections.push_back(x_limits);
    projections.axon_projections.push_back(y_limits);
    
}

void Axon::set_spheres(std::vector<Dynamic_Sphere> spheres_to_add, int axon_id){

    if (spheres_to_add.size() != 0){

        this->begin = spheres_to_add[0].center;

        this->spheres = spheres_to_add;

        this->end = spheres_to_add[spheres_to_add.size()-1].center;
    }
    // create projections
    add_projection(axon_id);

}

bool Axon::isNearAxon(Vector3d &position,  double distance_to_be_inside){
    bool isnear = false;
    Vector2d x_limits = projections.axon_projections[0];
    Vector2d y_limits = projections.axon_projections[1];

    if (position[0] >= x_limits[0]-distance_to_be_inside && position[0] <= x_limits[1]+distance_to_be_inside){
        if (position[1] >= y_limits[0]-distance_to_be_inside && position[1] <= y_limits[1]+distance_to_be_inside){
            return true;
        }
    }
    return false;
}

bool Axon::isPosInsideAxon(Vector3d &position,  double distance_to_be_inside, bool swell){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other axons -> check with max_radius so there is room for swelling

    if(isNearAxon(position, distance_to_be_inside)){
        std::vector<std::vector<Projections::projection_pt>> coliding_projs;
        double rad;
        if (swell){
            rad = max_radius;
            coliding_projs = projections_max.find_collisions_all_axes(position, rad, id);
        }
        else{
            rad = radius;
            coliding_projs = projections.find_collisions_all_axes(position, rad, id);
        }
    
        if (coliding_projs.size() == 3){ 

            // for all coliding objects in x 
            for(unsigned j = 0; j < coliding_projs[0].size() ; j++){ 
                const Projections::projection_pt coliding_proj = coliding_projs[0][j];
                // if the same coliding objects are also in y and z but are not from same objects
                bool colliding_all_axes;
                if (swell){
                    colliding_all_axes = (projections_max.isProjInside(coliding_projs[1], coliding_proj) && projections_max.isProjInside(coliding_projs[2], coliding_proj));
                }
                else{
                    colliding_all_axes = (projections.isProjInside(coliding_projs[1], coliding_proj) && projections.isProjInside(coliding_projs[2], coliding_proj));
                }
                if (colliding_all_axes){
                    Dynamic_Sphere sphere_ = spheres[coliding_proj.sph_id];

                    if (sphere_.distSmallerThan(position, distance_to_be_inside + rad)){ 
                        
                        return true;
                        break;
                    }  
                } 
            }
        }
    }
    return false;

} 


std::vector<double> Axon::Distances_to_Spheres(Vector3d pos){
    std::vector<double> distances;
    distances.clear();
    for (unsigned i=0; i< spheres.size(); ++i){
        //if (spheres[i].center[0] == begin[0] && spheres[i].center[1] == begin[1]){ 
        Vector3d m = pos - spheres[i].center;
        double distance_to_sphere = m.norm() - spheres[i].radius;
        distances.push_back(distance_to_sphere);
        //} 

    }
    return distances;
}

std::vector<double> Axon::Distances_to_Centers(Vector3d pos){
    std::vector<double> distances;
    distances.clear();
    for (unsigned i=0; i< spheres.size(); ++i){
        //if (spheres[i].center[0] == begin[0] && spheres[i].center[1] == begin[1]){ 
        Vector3d m = pos - spheres[i].center;
        double distance_to_sphere = m.norm();
        distances.push_back(distance_to_sphere);
        //} 

    }
    return distances;
}




std::vector<double> Axon::Distances_to_Spheres(Walker &w){

    Vector3d O;
    w.getVoxelPosition(O);
    return Distances_to_Spheres(O);

}


std::vector<Dynamic_Sphere> Axon::closestSpheres(Vector3d pos){

    Vector3d O = pos;
    std::vector <double> distances;
    std::vector <Dynamic_Sphere> closest_spheres;
    closest_spheres.clear();
    distances.clear();

    //string message = "spheres size : " +std::to_string(spheres.size())+"\n";
    //SimErrno::info(message,cout);
    for (unsigned i=0; i< spheres.size(); ++i){
        Vector3d m = O - spheres[i].center;
        double distance_to_sphere = m.norm();
        distance_to_sphere = distance_to_sphere - radius;
        distances.push_back(distance_to_sphere);
    }

    std::vector <double> d = distances;
    std::vector <Dynamic_Sphere> s = spheres;

    for (unsigned i=0; i< spheres.size(); ++i){
        auto it = std::min_element(std::begin(d), std::end(d));
        unsigned ind = std::distance(std::begin(d), it);
        closest_spheres.push_back(s[ind]);
        //string message = "closest_spheres size : " +std::to_string(closest_spheres.size())+"\n";
        //SimErrno::info(message,cout);
        s.erase(s.begin() + ind);
        d.erase(d.begin() + ind);
    }

    return closest_spheres;
}

std::vector<Dynamic_Sphere> Axon::closestSpheres(Walker &w){
    Vector3d O;
    w.getVoxelPosition(O);
    return closestSpheres(O);

}

double Axon::minDistance(Walker &w){

    std::vector<double> distances;
    distances.clear();
    distances = Distances_to_Spheres(w);
    double min = *std::min_element(std::begin(distances), std::end(distances ));

    return min;

}

double Axon::minDistance(Vector3d pos){

    std::vector<double> distances = Distances_to_Spheres(pos);
    double min = *std::min_element(std::begin(distances), std::end(distances ));
    return min;

}

double Axon::minDistanceCenter(Vector3d pos){

    std::vector<double> distances = Distances_to_Centers(pos);
    double min = *std::min_element(std::begin(distances), std::end(distances ));
    return min;

}

bool Axon::intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere &s, Vector3d &step, double &step_length, Vector3d &pos, bool isintra, double &c){
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    Vector3d m = pos - s.center;
    double rad = s.radius;

    double a = 1;
    double b = (m.dot(step));
    c = m.dot(m) - rad*rad;

    if(b > EPS_VAL && c > EPS_VAL)
        return false;
    
    double discr = b*b - a*c;

    if (discr < 0 ){
        return false;
    }

    t1 = (-b + sqrt(discr))/(a);
    t2 = (-b - sqrt(discr))/(a);

    return true;

}

double getAverage(std::vector<double> const& v) {
    if (v.empty()) {
        return NAN;
    }
    double sum = 0;
    for (unsigned i=0; i< v.size(); ++i){
        sum += v[i];
    }
 
    return sum / v.size();
}

bool Axon::checkCollision(Walker &walker, Vector3d &step, double &step_lenght, Collision &colision)
{
    string message;
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d next_step = step*step_lenght+O;

    bool isintra;
    bool next_step_intra = isPosInsideAxon(next_step, 0, false);
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
    //cout << "walker is intra : " << isintra_ << endl;
    //cout << "walker location is intra : " << (walker.location == Walker::intra) << endl;
    //cout << "walker isintra : " << isintra << endl;
    // if is intra and so is next step -> no collision with border
    if(isintra && next_step_intra){
        colision.type = Collision::null;
        cout << "wuuut" << endl;
        return false;
    }

    // distances to intersections
    std::vector<double> dist_intersections;
    std::vector<double> cs;
    std::vector<int> sph_ids;
    sph_ids.clear();
    dist_intersections.clear();
    int sph_id;

    
    for (unsigned i=0 ; i< spheres.size(); ++i){
        // distances to collision
        double t1;
        double t2;
        double c;
        bool intersect = intersection_sphere_vector(t1, t2, spheres[i], step, step_lenght, O, isintra, c);
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

        Dynamic_Sphere colliding_sphere = spheres[sphere_ind];

        colision.type = Collision::hit;
        colision.rn = cs[index_];

        if (!isintra){

            if(colision.rn <-1e-10){
                colision.col_location = Collision::inside;
                walker.in_obj_index = -1;
            }
            else if(colision.rn >1e-10){
                colision.col_location = Collision::outside;
            }
            else{
                colision.col_location = Collision::unknown;
            }
        }
        else{
            colision.col_location = Collision::inside;
            walker.in_obj_index = -1;

        }
            

        cout << "c :" << colision.rn << endl;
            

        colision.t = fmin(dist_to_collision,step_lenght);
        colision.colision_point = walker.pos_v + colision.t*step;

        //cout << "walker position : (" <<  walker.pos_v[0] << ", " << walker.pos_v[1] << "," << walker.pos_v[2] << ")" << endl;    

        //cout << "collision position : (" <<  colision.colision_point[0] << ", " << colision.colision_point[1] << "," << colision.colision_point[2] << ")" << endl;    

        //cout << "collision with sphere :" << sphere_ind << endl;
            
        //Normal point
        Vector3d normal = (colision.colision_point- colliding_sphere.center).normalized();
        Vector3d temp_step = step;
        elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);
        colision.bounced_direction = temp_step.normalized();
        return true;
        
    }
    else{
        //cout<< "nothing" << endl;
        colision.type = Collision::null;
        return false;   
    }

}



