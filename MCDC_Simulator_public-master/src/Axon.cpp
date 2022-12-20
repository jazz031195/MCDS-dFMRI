#include "Axon.h"
#include "dynamic_sphere.h"
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
    active_state = ax.active_state;


}

std::vector<double> Axon::Distances_to_Spheres(Eigen::Vector3d pos){
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



std::vector<double> Axon::Distances_to_Spheres(Walker &w){

    Vector3d O;
    w.getVoxelPosition(O);
    return Distances_to_Spheres(O);

}


std::vector<Dynamic_Sphere> Axon::closestSpheres(Eigen::Vector3d pos){

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

double Axon::minDistance(Eigen::Vector3d pos){

    std::vector<double> distances = Distances_to_Spheres(pos);
    double min = *std::min_element(std::begin(distances), std::end(distances ));
    return min;

}

double Axon::intersection_sphere_vector(Dynamic_Sphere &s, Eigen::Vector3d &step, double& step_length, Eigen::Vector3d &pos, bool isintra){
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    Vector3d m = pos - s.center;
    double rad = s.radius;

    double a = 1;
    double b = (m.dot(step));
    double c = m.dot(m) - rad*rad;
    
    double discr = b*b - a*c;


    if (discr < 0 ){
        if (isintra){
            return 0;
        } 
        else{ 
            return 10;
        } 
    }

    double t1 = (-b + sqrt(discr))/a;
    double t2 = (-b - sqrt(discr))/a;

    if (isintra == false){ 
        if (t1 >= 0 && t2 >= t1){

            return t1;
        }
        else if (t2 >= 0 && t1 >= t2){
            return t2;
            
        } 
        else {
            return 10;
        }
    }
    else {
        
        if (t2 >= t1 && t2 >= 0){
            
            return t2;
        }
        else if (t1 >= t2 && t1 >= 0){

            return t1;
        } 
        else {
            return 0;
        }

    }  

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

bool Axon::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    string message;
    Vector3d O;
    walker.getVoxelPosition(O);
    Eigen::Vector3d next_step = step*step_lenght+O;
    
    bool collision_check = false;
    if (walker.initial_location == Walker::intra){
        colision.col_location = Collision::inside;

        bool isinside_axon = (isInside(O, barrier_tickness));


        if (isinside_axon){
            

            if (isInside(next_step, -EPS_VAL)){
                colision.type = Collision::null;
                return false;
            }

            bool isintra = true;
            

            std::vector<double> dist_intersections;
            dist_intersections.clear();

            double t;
            for (unsigned i=0 ; i< spheres.size(); ++i){
                //double dist = (spheres[i].center-O).norm();
                //dist -= spheres[i].radius;
                //if (dist > barrier_tickness+step_lenght){
                //    dist_intersections.push_back(0);
                //}
                //else{

                    // position of intersections
                t = intersection_sphere_vector(spheres[i], step,step_lenght, O, isintra);
                //next_step = t*step+O;
                //if (isInside(next_step, barrier_tickness)){
                dist_intersections.push_back(t);
                t = 0;
                //}

                //}

            }
            auto max_distance_int = std::max_element(std::begin(dist_intersections), std::end(dist_intersections));


            if (*max_distance_int != 0 ){
                unsigned index_ = std::distance(std::begin(dist_intersections), max_distance_int);
                //message = "walker id : "+to_string(walker.index)+ ", max_distance_int: " +to_string(*max_distance_int)+" \n";
                //SimErrno::info(message,cout);
                collision_check = spheres[index_].checkCollision_(walker, step, step_lenght, colision, isintra);

                //if(collision_check){
                //    message = "is inside axon  : "+to_string(isInside(colision.colision_point, barrier_tickness))+ "  \n";
                //    SimErrno::info(message,cout);
                //}
                //else{
                //    next_step = step_lenght*step+O;
                //    message = "is inside axon  : "+to_string(isInside(next_step, barrier_tickness))+ "  \n";
                //    SimErrno::info(message,cout);

                //}
                return collision_check;
            
            }
            else{
                colision.type = Collision::null;
                return false;
            }   


        }

        else {
            colision.type = Collision::null;
            return false;
        }
    }


    else if (walker.initial_location == Walker::extra) {

        if (!isInside(next_step, barrier_tickness)){
            colision.type = Collision::null;
            return false;
        }


        bool isintra = false;

        // extracelluar particles

        // distances to intersections
        std::vector<double> dist_intersections;
        dist_intersections.clear();

        //string message = "first index : " +std::to_string(first_ind)+", last_index : "+std::to_string(last_ind)+" \n";
        //SimErrno::info(message,cout);
        double t;
        for (unsigned i=0 ; i< spheres.size(); ++i){

            //double dist = (spheres[i].center-O).norm();
            //dist -= spheres[i].radius;
            //if (dist > barrier_tickness+step_lenght){
                // spheres that are far away
                //dist_intersections.push_back(10);
            //}
            //else {
                // position of intersections
            t = intersection_sphere_vector(spheres[i], step,step_lenght, O, isintra);
            dist_intersections.push_back(t);
            t = 10;

            //}

        }

        auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));

        if (*min_distance_int != 10){
            unsigned index_ = std::distance(std::begin(dist_intersections), min_distance_int);

            return spheres[index_].checkCollision_(walker, step, step_lenght, colision, isintra);

        }
        else{
            colision.type = Collision::null;
            return false;
        }   

    }

}

bool Axon::isInside(Eigen::Vector3d pos, double distance_to_be_inside){
    Vector3d m;
    for (unsigned i=0; i< spheres.size(); ++i){
        m = pos - spheres[i].center;
        double distance_to_sphere = m.norm() - spheres[i].radius;
        if (distance_to_sphere <= distance_to_be_inside){
            return true;
            break;
        }
    
    }
    return false;
}

