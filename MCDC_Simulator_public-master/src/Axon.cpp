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
    Vector3d smallest_pos;
    Vector3d largest_pos;

    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {
        double smallest_pos_ = 100;
        double largest_pos_ = 0;

        for (size_t i=0; i < spheres.size(); ++i){ 
            // TODO : ask jasmine why 
            double position1, position2;//, position1_, position2_;
            int sph_id = i;

            // projections_max
            // center + radius
            position1 = spheres[i].center[axis] + spheres[i].max_radius;
            Projections::projection_pt p1 {position1, axon_id, sph_id};
            if (position1 < smallest_pos_){
                smallest_pos_ = position1;
            }
            if (position1 > largest_pos_){
                largest_pos_= position1;
            }
            // center - radius
            position2 = spheres[i].center[axis] - spheres[i].max_radius;
            Projections::projection_pt p2 {position2, axon_id, sph_id};

            if (position2 < smallest_pos_){
                smallest_pos_ = position2;
            }
            if (position2 > largest_pos_){
                largest_pos_= position2;
            }
            projections_max.append_right_place(p1,  p2, axis);

            // // projections
            // // center + radius
            // position1_ = spheres[i].center[axis] + spheres[i].radius;
            Projections::projection_pt p1_ {position1, axon_id, sph_id};
            // // center - radius
            // position2_ = spheres[i].center[axis] - spheres[i].radius;
            Projections::projection_pt p2_ {position2, axon_id, sph_id};

            projections.append_right_place(p1_,  p2_, axis);
        }

        smallest_pos[axis] = smallest_pos_;
        largest_pos[axis] = largest_pos_;
        
    }

    Vector2d x_limits = {smallest_pos[0], largest_pos[0]};
    Vector2d y_limits = {smallest_pos[1], largest_pos[1]};
    Vector2d z_limits = {smallest_pos[2], largest_pos[2]};
    projections_max.axon_projections.push_back(x_limits);
    projections_max.axon_projections.push_back(y_limits);
    projections_max.axon_projections.push_back(z_limits);
    projections.axon_projections.push_back(x_limits);
    projections.axon_projections.push_back(y_limits);
    projections.axon_projections.push_back(z_limits);
    
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

    size_t count_isnear = 0;
    for (size_t axis=0 ; axis < 3 ; axis++)
    {
        Vector2d axis_limits = projections.axon_projections[axis];

        if ((position[axis] >= axis_limits[0] - distance_to_be_inside) && 
            (position[axis] <= axis_limits[1] + distance_to_be_inside))
        {
            ++count_isnear;
        }
    }
    if (count_isnear == 3)
        return true;
    

    return false;
}

bool Axon::isPosInsideAxon(Vector3d &position,  double distance_to_be_inside, bool swell_, std::vector<int> sphere_ids){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other axons -> check with max_radius so there is room for swelling
    std::vector<std::vector<Projections::projection_pt>> coliding_projs;
    bool colliding_all_axes;
    Dynamic_Sphere sphere_ ;
    double rad;
    double rad_;
    // if position is in box with axon inside
    if(isNearAxon(position, distance_to_be_inside)){
        if (swell_){
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
                if (swell_){
                    colliding_all_axes = (projections_max.isProjInside(coliding_projs[1], coliding_proj) && projections_max.isProjInside(coliding_projs[2], coliding_proj));
                }
                else{
                    colliding_all_axes = (projections.isProjInside(coliding_projs[1], coliding_proj) && projections.isProjInside(coliding_projs[2], coliding_proj));
                }
                if (colliding_all_axes){
                    sphere_ = spheres[coliding_proj.sph_id];
                    if (swell_){
                        rad_ = sphere_.max_radius;
                    }
                    else{
                        rad_ = sphere_.radius;
                    }
                    if (sphere_.distSmallerThan(position, distance_to_be_inside + rad_)){ 

                        sphere_ids.push_back(coliding_proj.sph_id);
                        
                    }  
                } 
            }
        }
        if (sphere_ids.size() > 0){
            return true;
        }
        else{
            return false;
        }

    }
    //if (minDistance(position) < distance_to_be_inside){

    //    cout << "not working " << endl;
    //    if (isNearAxon(position, distance_to_be_inside)){
    //        if(coliding_projs.size() == 3){
    //            if(colliding_all_axes){
    //                if(!sphere_.distSmallerThan(position, distance_to_be_inside + rad)){
    //                    cout <<"distSmallerThan doesn't work" << endl;
    //                }
    //            }
    //            else{
    //                cout << "not colliding_all_axes" << endl;
    //            }
    //        }
    //        else{
    //            cout << "coliding_projs.size() != 3" << endl;
    //        }
    //    }
    //    else{
    //        cout << "not near" << endl;
    //    }
    //    return true;
    //}
    //else{
    //}
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

int Axon::closest_sphere_dichotomy(Walker &walker, double &step_lenght){
    Vector3d O;
    walker.getVoxelPosition(O);
    int number_spheres = spheres.size();
    int i=0;
    int i_last = number_spheres-1;
    int half_way;
    double first_distance;
    double last_distance;
    while ((i_last-i)>1){
        half_way = int((i_last+i)/2);
        first_distance =(spheres[i].center-O).norm() - spheres[i].radius;
        last_distance =(spheres[i_last].center-O).norm() - spheres[i_last].radius;
        if(first_distance < last_distance){
            i_last = half_way;
        } 
        else{
            i = half_way;
        } 
        count += 1;
    }
    if (first_distance < last_distance){
        return i;
    }
    else{
        return i_last;
    }  
 
}

Vector3d find_nearest_point_on_skeleton(Vector3d collision_point, Vector3d sphere_center1, Vector3d sphere_center2){
    double distance_to_point = (collision_point-sphere_center1).dot(sphere_center2-sphere_center1)/(sphere_center1-sphere_center2).norm();
    Vector3d point = sphere_center1+ (sphere_center2-sphere_center1)*distance_to_point;
    return point;
}


bool Axon::checkCollision(Walker &walker, Vector3d &step, double &step_lenght, Collision &colision)
{
    string message;
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d next_step = step*step_lenght+O;

    bool isintra;
    std::vector<int> col_sphere_ids_;
    bool next_step_is_intra;
    if (walker.location == Walker::intra){
        isintra = true;
        next_step_is_intra = isPosInsideAxon(next_step, -barrier_tickness, false, col_sphere_ids_);
    }
    else if (walker.location == Walker::extra) {
        isintra = false;
    }
    else{
        if(walker.initial_location == Walker::intra){
            isintra = true;
            next_step_is_intra = isPosInsideAxon(next_step, -barrier_tickness, false, col_sphere_ids_);
        }
        else{
            isintra = false;
        }
    }

    if(isintra && next_step_is_intra){
        colision.type = Collision::null;
        return false;
    }


    // is near axon (inside box)
    if (!isNearAxon(O, step_lenght+barrier_tickness)){
        colision.type = Collision::null;
        //cout << "not near axon" << endl;
        return false;
    }

    // is inside or near a sphere
    //std::vector<int> col_sphere_ids;
    //bool isnearspheres = isPosInsideAxon(O, step_lenght+barrier_tickness*10, false, col_sphere_ids);
    //if(!isnearspheres){
    //    cout << "not near spheres" << endl;
    //    colision.type = Collision::null;
    //    return false;
    //}


    // distances to intersections
    std::vector<double> dist_intersections;
    std::vector<double> cs;
    std::vector<int> sph_ids;
    sph_ids.clear();
    dist_intersections.clear();
    int sph_id;

 
    for (unsigned j=0 ; j < spheres.size(); ++j){

        int i = j;
        //int i = col_sphere_ids[j];
        // distances to collision
        double t1;
        double t2;
        double c;
        bool intersect = intersection_sphere_vector(t1, t2, spheres[i], step, step_lenght, O, isintra, c);
        // if ((!isintra && c > 0) || (isintra && c < 0)){
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
        // }
    }

    if(dist_intersections.size() > 0){
        //cout << "collision" << endl;
        unsigned index_ ;

        if (!isintra){
            auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
            index_ = std::distance(std::begin(dist_intersections), min_distance_int);
        }
        else{
            auto max_distance_int = std::max_element(std::begin(dist_intersections), std::end(dist_intersections));
            index_ = std::distance(std::begin(dist_intersections), max_distance_int);
        }

        const auto sphere_ind = static_cast<size_t>(sph_ids[index_]);

        double dist_to_collision = dist_intersections[index_];

        Dynamic_Sphere colliding_sphere = spheres[sphere_ind];

        colision.type = Collision::hit;
        colision.rn = cs[index_];

        if(!isintra)
        {
            if(colision.rn < -1e-10){
            colision.col_location = Collision::inside;
            }
            else if(colision.rn > 1e-10){
                colision.col_location = Collision::outside;
            }
            else{
                colision.col_location = Collision::unknown;
            }
        }
        else if (isintra && !check_negatif(cs))
        {
            colision.col_location = Collision::outside; 
        }
        

        colision.t = fmin(dist_to_collision,step_lenght);
        colision.colision_point = walker.pos_v + colision.t*step;

        //Normal point
        Vector3d normal;
        const size_t left_index  = sphere_ind > 0 ? sphere_ind - 1 : 0;
        const size_t right_index = sphere_ind + 1 < spheres.size() ? sphere_ind + 1 : spheres.size() - 1;
        const double d1 = (colision.colision_point-spheres[left_index].center).norm();
        const double d2 = (colision.colision_point-spheres[right_index].center).norm();
        double d;
        Vector3d close_sphere_pos;
        if (d1 < d2){
            close_sphere_pos = spheres[left_index].center;
            d = d1;
        }
        else{
            close_sphere_pos = spheres[right_index].center;
            d = d2;
        }

        // walker is at intersection of spheres
        if (abs((colision.colision_point - spheres[sphere_ind].center).norm() - d) < EPS_VAL){
            Vector3d point_on_skeleton = find_nearest_point_on_skeleton(colision.colision_point, spheres[sphere_ind].center, close_sphere_pos);
            normal = (colision.colision_point - point_on_skeleton).normalized();
        }
        else{
            normal = (colision.colision_point - spheres[sphere_ind].center).normalized();
        }


        Vector3d temp_step = step;
        elasticBounceAgainsPlane(walker.pos_v, normal, colision.t, temp_step);
        colision.bounced_direction = temp_step.normalized();

        return true;
        
    }
    else{
        //cout << "no collision" << endl;
        colision.type = Collision::null;
        return false;   
    }

}

bool Axon::check_negatif(vector<double> list)
{
    for (size_t i=0 ; i < list.size(); ++i)
    {
        if (list[i] < 1e-10)
            return true;
    }
    return false;
}

double Axon::volumeAxon()
{
    double volume = 0;
    // double tortuosity;
    double ax_length = 0;
    double mean_rad  = 0;
    if (spheres.size() > 0)
    {
        for (uint j = 0; j < spheres.size(); j++){
            if (j > 0){
                // Length between two adjacent spheres' centers
                double l = (spheres[j-1].center - spheres[j].center).norm();
                ax_length += l;
            }
            mean_rad += spheres[j].radius;
        }
        mean_rad   = mean_rad/spheres.size();
        tortuosity = ax_length/((this->begin - this->end).norm());
        volume     = M_PI * mean_rad * mean_rad * ax_length;
    }
    else
    {
        //TODO : throw an error
    }
    
    return volume;
}


