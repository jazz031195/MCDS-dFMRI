#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include "simerrno.h"

using namespace Eigen;
using namespace std;


Axon::Axon(const Axon &ax)
{

    swell = ax.swell;
    id = ax.id;
    spheres = ax.spheres;
    radius = ax.radius;
    begin = ax.begin;
    end = ax.end;
    volume_inc_perc = ax.volume_inc_perc;
    projections = ax.projections;
    min_radius = ax.min_radius;

}

void Axon::add_projection(){

    // 2d limits of axon
    Vector2d x_limits;
    Vector2d y_limits;
    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {

        for (int i=0; i<spheres.size(); ++i){ 

            double position1, position2;

            // projections
            // center + radius
            position1 = spheres[i].center[axis] + spheres[i].radius;
            Projections::projection_pt p1 {position1, id, i};
            // center - radius
            position2 = spheres[i].center[axis] - spheres[i].radius;
            Projections::projection_pt p2 {position2, id, i};

            projections.append_right_place(p1,  p2, axis);
        }

        if (axis == 0 && projections.sph_projections_x.size()>0) {
            x_limits = {projections.sph_projections_x[0].position, projections.sph_projections_x[projections.sph_projections_x.size()-1].position};
        }
        else if (axis == 1 && projections.sph_projections_y.size()>0) {
            y_limits = {projections.sph_projections_y[0].position, projections.sph_projections_y[projections.sph_projections_y.size()-1].position};
        }
    }
    if (x_limits.size()>0){
        projections.axon_projections.push_back(x_limits);
    }
    if (y_limits.size()>0){
        projections.axon_projections.push_back(y_limits);
    }

}

void Axon::set_spheres(std::vector<Dynamic_Sphere> spheres_to_add){

    // set axon_id of spheres to the id of axon
    for (int i=0; i<spheres_to_add.size(); ++i){
        spheres_to_add[i].ax_id = id;
        spheres_to_add[i].id = i;
    }
    if (spheres_to_add.size() != 0){

        this->begin = spheres_to_add[0].center;

        this->spheres = spheres_to_add;

        this->end = spheres_to_add[spheres_to_add.size()-1].center;

        // create projections
        add_projection();
    }



}

bool check_with_edge(Vector3d position, Vector2d x_limits, Vector2d y_limits){ 
    if ((position[0] >=  x_limits[0])  && (position[0] <= x_limits[1])){
        if ((position[1] >= y_limits[0]) && (position[1] <=  y_limits[1])){
            return true;
        }
    }
    return false;
} 

bool Axon::isNearAxon(Vector3d &position,  double distance_to_be_inside){
    bool isnear = false;
    Vector2d x_limits = projections.axon_projections[0];
    Vector2d y_limits = projections.axon_projections[1];

    if (distance_to_be_inside < 0){
        distance_to_be_inside = 0;
    } 
    x_limits ={x_limits[0]-distance_to_be_inside , x_limits[1]+distance_to_be_inside };
    y_limits ={y_limits[0]-distance_to_be_inside , y_limits[1]+distance_to_be_inside } ;

    if(check_with_edge(position, x_limits, y_limits)){
        return true;
    }  

    if (x_limits[0]<0){
        if(check_with_edge(position,{x_limits[0]+ end[2], x_limits[1] + end[2]} , y_limits)){
            return true;
        } 
    }
    else if (x_limits[0]>end[2]){
        if(check_with_edge(position,{x_limits[0] - end[2], x_limits[1] - end[2]} , y_limits)){
            return true;
        } 
    } 

    if (y_limits[0]<0){
        if(check_with_edge(position,x_limits,{y_limits[0]+ end[2], y_limits[1] + end[2]} )){
            return true;
        } 
    }
    else if (y_limits[0]>end[2]){
        if(check_with_edge(position,x_limits, {y_limits[0] - end[2], y_limits[1] - end[2]})){
            return true;
        } 
    } 

    return false;
}

bool Axon::isPosInsideAxon(Vector3d &position,  double distance_to_be_inside, std::vector<int>& sphere_ids){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other axons -> check with max_radius so there is room for swelling
    std::vector<std::vector<Projections::projection_pt>> coliding_projs;
    bool colliding_all_axes;
    Dynamic_Sphere sphere_ ;
    double rad;
    sphere_ids.clear();

    // if position is in box with axon inside
    if(isNearAxon(position, distance_to_be_inside)){
        // find all projections in between the two projections of the edges

        rad = radius+1000*barrier_tickness + distance_to_be_inside;
        
        coliding_projs = projections.find_collisions_all_axes(position, rad, id, distance_to_be_inside);

        if (coliding_projs.size() == 3){ 

            // for all coliding objects in x 
            for(unsigned j = 0; j < coliding_projs[0].size() ; j++){ 

                const Projections::projection_pt coliding_proj = coliding_projs[0][j];
                // if the same coliding objects are also in y and z but are not from same objects

                colliding_all_axes = (projections.isProjInside(coliding_projs[1], coliding_proj) && projections.isProjInside(coliding_projs[2], coliding_proj));

                if (colliding_all_axes){
                    sphere_ = spheres[coliding_proj.sph_id];
                    
                    if (sphere_.minDistance(position) < distance_to_be_inside){ 

                        sphere_ids.push_back(coliding_proj.sph_id);
                        
                        /*
                        cout << " Axon : "<< id <<"Position : [" << position[0] << ", "<< position[1] << ", " << position[2] << "]" << endl; 
                        cout << "           distance to sphere :" << coliding_proj.sph_id << " : " <<  sphere_.minDistance(position) ;
                        cout << ", Sphere position : [" << sphere_.center[0] << ", "<< sphere_.center[1] << ", " << sphere_.center[2] << "]" << endl;  
                        */
                    }  
                } 
            }
        }
        if (sphere_ids.size() > 0){
            // sort by value
            sort(sphere_ids.begin(), sphere_ids.end());
            return true;
        }
        else{
            //cout << " Not inside axon :" << id << endl;
            return false;
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

bool Axon::intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere &s, Vector3d &step, double &step_length, Vector3d &pos, double &c){
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



Vector3d find_nearest_point_on_skeleton(Vector3d collision_point, Vector3d sphere_center1, Vector3d sphere_center2){
    double distance_to_point = (collision_point-sphere_center1).dot(sphere_center2-sphere_center1)/(sphere_center1-sphere_center2).norm();
    Vector3d point = sphere_center1+ (sphere_center2-sphere_center1)*distance_to_point;
    return point;
}

bool check_inside(vector<double> list_c){
    for (unsigned i=0 ; i< list_c.size(); ++i){
        if (list_c[i] < 1e-10){
            return true;
            break;        
        }
    }
    return false;
}

bool check_outside(vector<double> list_c){
    for (unsigned i=0 ; i< list_c.size(); ++i){
        if (list_c[i] < -1e-10){
            return false;      
        }
    }
    return true;
}

bool check_on_edge(vector<double> list_c){
    for (unsigned i=0 ; i< list_c.size(); ++i){
        if (list_c[i] < 1e-10 && list_c[i] > -1e-10){
            return true;
            break;        
        }
    }
    return false;
}

bool Axon::checkCollision(Walker &walker,  Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    
    string message;
    Eigen::Vector3d O;
    walker.getVoxelPosition(O);

    // distances to intersections
    std::vector<double> dist_intersections;
    // values indicating whether the walker is inside or outside a sphere
    std::vector<double> cs;

    std::vector<int> sph_ids;
    std::vector<double> all_cs;

    int first_ind;
    int last_ind;

    first_ind = 0;
    last_ind = spheres.size();
    

    
    for (unsigned i= first_ind ; i< last_ind; ++i){
        
        // distances to collision
        double t1;
        double t2;
        double c;
        bool intersect = intersection_sphere_vector(t1, t2, spheres[i], step, step_lenght, O, c); 
        int non_intersecting = 0;
        double limit_length = -0.5;
        bool has_intersected= false;

        if (intersect){
            non_intersecting = 0;
            has_intersected= true;

            all_cs.push_back(c);


            bool condition;

            condition = true;

            if (condition){ 

                //if the collision are too close or negative.
                if(Walker::bouncing){
                    if( t1 >= EPS_VAL && t1 <= step_lenght + barrier_tickness){

                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
                else{
                    if( t1 >= 0 && t1 <= step_lenght + barrier_tickness){

                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }  
            }
            condition = true;
            //condition = true;
            //if the new position is at edge of axon, at the edge of sphere[i] but inside the neighbours
            if  (condition){ 

                if (Walker::bouncing){
                    if(t2 >= EPS_VAL && t2 <= step_lenght + barrier_tickness){

                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
                else{
                    if(t2 >= 0 && t2 <= step_lenght + barrier_tickness){

                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }

            }
            
        } 
        else{
            non_intersecting +=1;
        }
        if(has_intersected and non_intersecting > 5){
            break;
        }

    }

    if(dist_intersections.size() > 0){

        auto min_distance_int = std::max_element(std::begin(dist_intersections), std::end(dist_intersections));
        if (walker.initial_location== Walker::intra){
            min_distance_int = std::max_element(std::begin(dist_intersections), std::end(dist_intersections));
        }
        else{
            min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
        }
        unsigned index_ = std::distance(std::begin(dist_intersections), min_distance_int);
        
        int sphere_ind = sph_ids[index_];

        double dist_to_collision = dist_intersections[index_];

        colision.type = Collision::hit;
        colision.rn = cs[index_];  
        
 
        std::cout << "Collision, sphere index :" << sph_ids[index_] << endl;
        std::cout << "Collision, walker position :" << O << endl;
        std::cout << "Collision, collision position :" << O+step*dist_to_collision << endl;
        
        colision.t = fmin(dist_to_collision,step_lenght);
        colision.colision_point = walker.pos_v + colision.t*step;
        Eigen::Vector3d normal = (colision.colision_point- spheres[sphere_ind].center).normalized();
        
        Eigen::Vector3d temp_step = step;
        if (walker.initial_location== Walker::extra){
            colision.col_location = Collision::outside;

            elasticBounceAgainsPlane_extra(walker.pos_v,normal,colision.t,temp_step);
        }
        else{
            colision.col_location = Collision::inside;

            elasticBounceAgainsPlane_intra(walker.pos_v,normal,colision.t,temp_step);
        }
        colision.bounced_direction = temp_step.normalized();
        colision.collision_objects = {sphere_ind, id};

        return true;
        
    }
    else{

        colision.type = Collision::null;
        return false;   
    }

}





