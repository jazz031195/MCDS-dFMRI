#include "dynamic_sphere.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include "simerrno.h"
using namespace Eigen;
using namespace std;



Dynamic_Sphere::Dynamic_Sphere(const Dynamic_Sphere &sph)
{
    center = sph.center;
    radius = sph.radius;
    swell = sph.swell;
    volume_inc_perc = sph.volume_inc_perc; 
    ax_id = sph.ax_id;
    active_state = sph.active_state;
    min_radius = sph.min_radius;
    max_radius = sph.max_radius;
    id = sph.id;

}

void Dynamic_Sphere::set_center(Eigen::Vector3d center_)
{
    this->center = center_; 

}

bool Dynamic_Sphere::checkCollision_(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision, bool &isintra)
{
    //tmax = step_length

    //Origin of the ray
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d m = O - this->center;

    // total distance
    double distance_to_sphere = m.norm();
    // collision distance
    double d_ = distance_to_sphere - max_radius;

    //If the minimum distance from the walker to the cylinder is more than
    // the actual step size, we can discard this collision.
    if(d_> EPS_VAL){
        if(d_ > step_lenght+barrier_tickness){
            colision.type = Collision::null;
            return false;
        }
    }

    double a = 1;
    double b = m.dot(step);
    double c = m.dot(m) - max_radius*max_radius;
    if (isintra == false){
        if(b > EPS_VAL && c > EPS_VAL){
            colision.type = Collision::null;
            return false;
        }
    }


    double discr = b*b - a*c;
    // no real roots
    if(discr <= 0.0){
        colision.type = Collision::null;
        return false;
    }

    //if we arrived here we need to compute the quadratic equation.
    return handleCollition(walker,colision,step,a,b,c,discr,step_lenght, isintra);

}


inline bool Dynamic_Sphere::handleCollition(Walker& walker, Collision &colision, Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length, bool& isintra){

    double t1 = (-b - sqrt(discr))/a;

    double t2 = (-b + sqrt(discr))/a;

    // take the longest path if walker is in axon
    if (isintra){
        colision.col_location = Collision::inside;
        if (t1> t2){
            t2 = t1;
        }
        else {
            t1 = t2;
        }
    }


    //if we are completely sure that no collision happened
    if( ( (t1 < 0.0) || (t1 > step_length+barrier_tickness) ) && ( (t2 < 0.0) || (t2 > step_length+barrier_tickness)) ){
        colision.type = Collision::null;
        return false;
    }

    // a spin that's bouncing ignores collision at 0 (is in a wall)
    if(walker.status == Walker::bouncing){
        //string message = "  Bouncing : axon id :"+to_string(ax_id)+", pos of sphere (" +std::to_string(center[0])+", "+std::to_string(center[1])+", "+std::to_string(center[2])+") \n";
        //SimErrno::info(message,cout);

        //if the collision are too close or negative.
        if( ( (t1 < EPS_VAL) || (t1 > step_length+barrier_tickness)) && (( t2 < EPS_VAL) || (t2 > step_length+barrier_tickness)) ){
            //message = "t too close or negative \n";
            //SimErrno::info(message,cout);
            colision.type = Collision::null;
            return false;
        }

        if( t1 >= EPS_VAL && t1 < t2)
            colision.t = fmin(t1,step_length);
        else
            colision.t = fmin(t2,step_length);
    }
    else{
        if( t1>0.0 && t1 <t2)
            colision.t = fmin(t1,step_length);
        else
            colision.t = fmin(t2,step_length);
    }


    colision.type = Collision::hit;
    colision.obstacle_ind = -1;

    if (isintra == false){

        if(c<-1e-10){
            colision.col_location = Collision::inside;

            //string message = "Inside ! : axon id :"+to_string(ax_id)+", pos of sphere (" +std::to_string(center[0])+", "+std::to_string(center[1])+", "+std::to_string(center[2])+") \n";
            //SimErrno::error(message,cout);
        //    walker.in_obj_index = -1;
        }
    else if(c>1e-10){

            colision.col_location = Collision::outside;
            
        }
        else{
            colision.col_location = Collision::unknown;
        }
    }

    colision.rn = c;
    colision.colision_point = walker.pos_v + colision.t*step;

    //Normal point
    Eigen::Vector3d normal = (colision.colision_point-this->center).normalized();
    Eigen::Vector3d temp_step = step;
    elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);
    colision.bounced_direction = temp_step.normalized();

    //string message = "Colision at t : "+to_string(colision.t)+" \n";
    //SimErrno::error(message,cout);    

    return true;

}

double Dynamic_Sphere::minDistance(Walker &w){

    //Origin of the ray
    Vector3d O;
    w.getVoxelPosition(O);
    Vector3d m = O - this->center;
    // minimum distance to the cylinder axis.
    double distance_to_sphere = m.norm();

    //Minimum distance to the shpere wall.
    double d_ = (distance_to_sphere - radius);
   // return d_>0.0?d_:0.0;
    return d_;
}

bool Dynamic_Sphere::isInside(Walker &w){

    //Minimum distance to the sphere wall.
    double d_ = minDistance(w);
   // return d_>0.0?d_:0.0;
    return d_ <= 0;
}

bool Dynamic_Sphere::isInside(Eigen::Vector3d pos, double distance_to_be_inside){
    double d_ = (pos - this->center).norm();
 
    d_ = d_-this->radius;
    
   // return d_>0.0?d_:0.0;
    return d_ <= distance_to_be_inside;
}

bool Dynamic_Sphere::distSmallerThan(Eigen::Vector3d pos, double distance){
    double d_ = (pos - this->center).norm();
    
    return d_ <= distance;
}
