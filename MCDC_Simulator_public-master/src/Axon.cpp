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
    distal_branching = ax.distal_branching;
    proximal_branching = ax.proximal_branching;
}

void Axon::add_projection(){

    // 2d limits of axon
    Vector2d x_limits;
    Vector2d y_limits;
    Vector2d z_limits;
    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {

        for (int i=0; i < static_cast<int>(spheres.size()); ++i){ 

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

        if (axis == 0) {
            x_limits = {projections.sph_projections_x[0].position, projections.sph_projections_x[projections.sph_projections_x.size()-1].position};
        }
        else if (axis == 1) {
            y_limits = {projections.sph_projections_y[0].position, projections.sph_projections_y[projections.sph_projections_y.size()-1].position};
        }
        else if (axis == 2) {
            z_limits = {projections.sph_projections_z[0].position, projections.sph_projections_z[projections.sph_projections_z.size()-1].position};
        }
    }

    projections.axon_projections.push_back(x_limits);
    projections.axon_projections.push_back(y_limits);
    projections.axon_projections.push_back(z_limits);

    /*
    cout << "Creating projections : " << endl;
    cout <<"x: "<<endl;
    for (int i=0; i<5; ++i){ 
        cout << projections.sph_projections_x[i].position << ", ";
    }
    cout << endl;
    cout <<"y: "<<endl;
    for (int i=0; i<5; ++i){ 
        cout << projections.sph_projections_y[i].position << ", ";
    }
    cout << endl;
    cout <<"z: "<<endl;
    for (int i=0; i<5; ++i){ 
        cout << projections.sph_projections_z[i].position << ", ";
    }
    cout << endl;
    */

}

void Axon::set_spheres(std::vector<Dynamic_Sphere>& spheres_to_add){

    // set axon_id of spheres to the id of axon
    for (size_t i=0; i<spheres_to_add.size(); ++i){
        spheres_to_add[i].ax_id = id;
        spheres_to_add[i].id = i;
    }
    if (spheres_to_add.size() != 0){

        this->begin = spheres_to_add[0].center;

        this->spheres = spheres_to_add;

        this->end = spheres_to_add[spheres_to_add.size()-1].center;
    }

    // create projections
    add_projection();

}

bool check_with_edge(Vector3d const& position, Vector2d const& x_limits, Vector2d const& y_limits)  
{ 
    if ((position[0] >=  x_limits[0])  && (position[0] <= x_limits[1])){
        if ((position[1] >= y_limits[0]) && (position[1] <=  y_limits[1])){
            return true;
        }
    }
    return false;
} 

bool check_with_edge(Vector3d const& position, Vector2d const& x_limits, Vector2d const& y_limits, Vector2d const& z_limits) 
{   
    if ((position[0] >= x_limits[0]) && (position[0] <= x_limits[1]) &&
        (position[1] >= y_limits[0]) && (position[1] <= y_limits[1]) &&
        (position[2] >= y_limits[0]) && (position[2] <= y_limits[1]))
          return true;
        
    return false;
} 

bool Axon::isNearAxon(Vector3d const&position,  double const& distance_to_be_inside){
    // bool isnear = false;
    Vector2d x_limits = projections.axon_projections[0];
    Vector2d y_limits = projections.axon_projections[1];
    Vector2d z_limits = projections.axon_projections[2];

    double dist;
    if (distance_to_be_inside < 0)
        dist = 0;
    else
        dist = distance_to_be_inside;

    x_limits ={x_limits[0]-dist , x_limits[1]+dist };
    y_limits ={y_limits[0]-dist , y_limits[1]+dist } ;
    z_limits ={z_limits[0]-dist , z_limits[1]+dist } ;

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

bool Axon::isPosInsideAxon(Vector3d const& position,  double const& distance_to_be_inside, vector<int>& sphere_ids, vector<double>& distances){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other axons -> check with max_radius so there is room for swelling
    std::vector<std::vector<Projections::projection_pt>> coliding_projs;
    bool colliding_all_axes;
    Dynamic_Sphere sphere_ ;
    double rad;
    sphere_ids.clear();
    distances.clear();

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
                        distances.push_back((sphere_.center - position).norm());
                    }  
                } 
            }
        }
        if (sphere_ids.size() > 0){
            double min_dist_to_center = 1000;
            int min_id = 1000;
            for(size_t sph_id=0; sph_id < sphere_ids.size(); ++ sph_id)
            {
                if (distances[sph_id] < min_dist_to_center)
                {
                    min_dist_to_center = distances[sph_id];
                    min_id = sphere_ids[sph_id];
                }
                    
            }
            sphere_ids.clear();
            sphere_ids.push_back(min_id);
            return true;
        }
        else
            return false;
    }

    return false;
} 

bool Axon::isPosInsideAxon(Vector3d const& position,  double const& distance_to_be_inside, vector<int>& sphere_ids){
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
                        
                        
                        // cout << " Axon : "<< id <<"Position : [" << position[0] << ", "<< position[1] << ", " << position[2] << "]" << endl; 
                        // cout << "           distance to sphere :" << coliding_proj.sph_id << " : " <<  sphere_.minDistance(position) ;
                        // cout << ", Sphere position : [" << sphere_.center[0] << ", "<< sphere_.center[1] << ", " << sphere_.center[2] << "]" << endl;  
                        
                    }  
                } 
            }
        }
        if (sphere_ids.size() > 0){
            // sort by value
            sort(sphere_ids.begin(), sphere_ids.end());
            double min_dist_to_center = 1000;
            int min_id = 1000;
            for(size_t sph_id=0; sph_id < sphere_ids.size(); ++ sph_id)
            {
                double dist_tmp = (spheres[sphere_ids[sph_id]].center - position).norm();
                cout << "dist " << sphere_ids[sph_id] << " " << dist_tmp << " " << sphere_ids[sph_id] << endl;
                if (dist_tmp < min_dist_to_center)
                {
                    min_dist_to_center = dist_tmp;
                    min_id = sphere_ids[sph_id];
                }
                    
            }
            cout << "min dist" << min_dist_to_center << endl;
            sphere_ids.clear();
            sphere_ids.push_back(min_id);
            // cout << "sph_ids size" << sphere_ids.size() << endl;
            // for(size_t i=0; i < spheres[sphere_ids[0]].neighboring_spheres.size(); ++i)
            //     cout << (spheres[sphere_ids[0]].neighboring_spheres[i]->center - position).norm() << endl;
            return true;
        }
        else{
            //cout << " Not inside axon :" << id << endl;
            return false;
        }
        

    }

    return false;
} 

std::vector<double> Axon::Distances_to_Spheres(Vector3d const& pos) const {
    std::vector<double> distances;
    distances.clear();
    for (unsigned i=0; i < spheres.size(); ++i){
        //if (spheres[i].center[0] == begin[0] && spheres[i].center[1] == begin[1]){ 
        Vector3d m = pos - spheres[i].center;
        double distance_to_sphere = m.norm() - spheres[i].radius;
        distances.push_back(distance_to_sphere);
        //} 

    }
    return distances;
}

std::vector<double> Axon::Distances_to_Centers(Vector3d const& pos) const {
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

double Axon::volumeAxon() const
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
        // tortuosity = ax_length/((this->begin - this->end).norm());
        volume     = M_PI * mean_rad * mean_rad * ax_length;
    }
    else
    {
        //TODO [ines]: throw an error
    }
    
    return volume;
}


std::vector<double> Axon::Distances_to_Spheres(Walker const& w) const{

    Vector3d O;
    w.getVoxelPosition(O);
    return Distances_to_Spheres(O);

}


std::vector<Dynamic_Sphere> Axon::closestSpheres(Vector3d const& pos) const{

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

std::vector<Dynamic_Sphere> Axon::closestSpheres(Walker const& w) const{
    Vector3d O;
    w.getVoxelPosition(O);
    return closestSpheres(O);

}

double Axon::minDistance(Walker const& w) const{

    std::vector<double> distances;
    distances.clear();
    distances = Distances_to_Spheres(w);
    double min = *std::min_element(std::begin(distances), std::end(distances ));

    return min;

}

double Axon::minDistance(Vector3d const& pos) const{

    std::vector<double> distances = Distances_to_Spheres(pos);
    double min = *std::min_element(std::begin(distances), std::end(distances ));
    return min;

}

double Axon::minDistanceCenter(Vector3d const& pos) const{

    std::vector<double> distances = Distances_to_Centers(pos);
    double min = *std::min_element(std::begin(distances), std::end(distances ));
    return min;

}
bool Axon::intersection_sphere_vector(double &t1, double &t2, Dynamic_Sphere const&s, Vector3d const&step, double const&step_length, 
                                      Vector3d const&pos, double &c){
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



Vector3d find_nearest_point_on_skeleton(Vector3d const& collision_point, Vector3d const& sphere_center1, Vector3d const& sphere_center2){
    double distance_to_point = (collision_point - sphere_center1).dot(sphere_center2 - sphere_center1)/(sphere_center1 - sphere_center2).norm();
    Vector3d point = sphere_center1 + (sphere_center2 - sphere_center1)*distance_to_point;
    return point;
}

bool check_inside(vector<double> const& list_c){
    for (unsigned i=0 ; i < list_c.size(); ++i){
        if (list_c[i] < 1e-10) 
            return true;
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


bool Axon::checkCollision(Walker &walker,  Vector3d const& step, double const& step_lenght, Collision &colision) 
{
    
    string message;
    Vector3d O;
    walker.getVoxelPosition(O);

    // distances to intersections
    std::vector<double> dist_intersections;
    // values indicating whether the walker is inside or outside a sphere
    std::vector<double> cs;

    std::vector<int> sph_ids;
    std::vector<double> all_cs;

    size_t first_ind = walker.in_sph_index[0]-10;
    if (first_ind <0) {
        first_ind = 0;
    }
    size_t last_ind = walker.in_sph_index[walker.in_sph_index.size()-1]+10;
    if (last_ind > spheres.size()) {
        last_ind = spheres.size();
    }
    /*
    cout << " ----- " << endl;
    
    for (unsigned k= 0 ; k< walker.in_sph_index.size(); ++k){
        cout << "   walker.in_sph_index [" << k << "] :" << walker.in_sph_index[k] << endl;
    }
    cout << "first_ind : " << first_ind <<", last_ind :" << last_ind << endl;
    */
    
    for (size_t i= first_ind ; i< last_ind; ++i){
        

        //cout << "Sphere : " << i << ", sph_id :" << spheres[i].id <<", Ax :" << id << endl;

        // distances to collision
        double t1;
        double t2;
        double c;
        bool intersect = intersection_sphere_vector(t1, t2, spheres[i], step, step_lenght,
                                                    O, c); 

        double limit_length = -EPS_VAL;

        if (intersect){

            all_cs.push_back(c);
            // check if the walker has previsouly collided with this sphere
            // bool same_colliding_object = false;
            std::vector<int> last_collision;
            walker.getLastCollision(last_collision);


            //if ((last_collision).size() == 2) {
                //cout << "Last collision of walker : sph : " << walker.last_collision[0] << ", ax : " << walker.last_collision[1] << endl;
            
            //    if ((i == last_collision[0])&&(id == last_collision[1])) {
            //        same_colliding_object = true;
                    //cout << "Same collision as previous !" << endl;
            //    }
            //}
            //cout << "Intersects with sphere " << j <<", t1 :" << t1 << ", t2 : " << t2 <<", c :" << c << endl;
            //if the new position is at edge of axon, at the edge of sphere[i] but inside the neighbours 
            bool condition;
            if (i==0){
                condition = !(spheres[i+1].isInside(step*t1+O, limit_length));
            }  
            else if ( i == spheres.size()){
                condition = !(spheres[i-1].isInside(step*t1+O, limit_length));
            } 
            else{
                condition = !(spheres[i+1].isInside(step*t1+O, limit_length) || spheres[i-1].isInside(step*t1+O, limit_length));
                //cout <<  "spheres[" << i+1 << "].minDistance(step*t1+O) : "<<spheres[i+1].minDistance(step*t1+O) << endl;
                //cout <<  "spheres[" << i-1 << "].minDistance(step*t1+O) : "<<spheres[i-1].minDistance(step*t1+O) << endl;
                //cout << "this sphere : spheres[" << i << "].minDistance(step*t1+O) : "<<spheres[i].minDistance(step*t1+O) << endl;
                //cout << "condition : " << condition << endl;
            }
            //condition = true; 
            if (condition){ 
                        
                //if the collision are too close or negative.
                if(Walker::bouncing){
                    if( t1 >= EPS_VAL && t1 <= step_lenght + barrier_tickness){
                        //cout << "   intersection, sphere (bouncing):" << i << endl;
                        //cout << "       c :" << c << endl;
                        //cout << "       t :" << t1 << endl;
                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
                else{
                    if( t1 >= 0 && t1 <= step_lenght + barrier_tickness){
                        //cout << "   intersection, sphere (bouncing):" << i << endl;
                        //cout << "       c :" << c << endl;
                        //cout << "       t :" << t1 << endl;
                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }  
            }
            if (i==0){
                condition = !(spheres[i+1].isInside(step*t2+O, limit_length));
            }  
            else if (i == spheres.size()){
                condition = !(spheres[i-1].isInside(step*t2+O, limit_length));
            } 
            else{
                condition = !(spheres[i+1].isInside(step*t2+O, limit_length) || spheres[i-1].isInside(step*t2+O, limit_length));
                //cout <<  "spheres[" << i+1 << "].minDistance(step*t1+O) : "<<spheres[i+1].minDistance(step*t2+O) << endl;
                //cout <<  "spheres[" << i-1 << "].minDistance(step*t1+O) : "<<spheres[i-1].minDistance(step*t2+O) << endl;
                //cout << "this sphere : spheres[" << i << "].minDistance(step*t1+O) : "<<spheres[i].minDistance(step*t2+O) << endl;
                //cout << "condition : " << condition << endl;
            
            }
            //condition = true; 
            //if the new position is at edge of axon, at the edge of sphere[i] but inside the neighbours
            if  (condition){ 
                if (Walker::bouncing){
                    if(t2 >= EPS_VAL && t2 <= step_lenght + barrier_tickness){
                        //cout << "   intersection, sphere :" << i << endl;
                        //cout << "       c :" << c<< endl;
                        //cout << "       t :" << t2 << endl;
                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }
                else{
                    if(t2 >= 0 && t2 <= step_lenght + barrier_tickness){
                        //cout << "   intersection, sphere :" << i << endl;
                        //cout << "       c :" << c<< endl;
                        //cout << "       t :" << t2 << endl;
                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                    }
                }

            }
            
        } 

    }

    if(dist_intersections.size() > 0){
        //cout << "dist_intersections.size() :" << dist_intersections.size() << endl;

        auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
        unsigned index_ = std::distance(std::begin(dist_intersections), min_distance_int);
        
        int sphere_ind = sph_ids[index_];

        double dist_to_collision = dist_intersections[index_];

        colision.type = Collision::hit;
        colision.rn = cs[index_];

        if (walker.initial_location== Walker::intra){ 
            if (check_inside(all_cs)){
                colision.col_location = Collision::inside;
            }
            else{
                colision.col_location = Collision::outside;
            } 
        } 
        else if(walker.initial_location== Walker::extra){
            if (check_outside(all_cs)){
                colision.col_location = Collision::outside;
            }
            else{
                colision.col_location = Collision::inside;
            } 

        }  
        
        //cout << "Collision, sphere :" << sph_ids[index_] << endl;
        colision.t = fmin(dist_to_collision,step_lenght);
        colision.colision_point = walker.pos_v + colision.t*step;
        Vector3d normal = (colision.colision_point- spheres[sphere_ind].center).normalized();
        

        Vector3d temp_step = step;
        elasticBounceAgainsPlane_intra(walker.pos_v,normal,colision.t,temp_step);
        colision.bounced_direction = temp_step.normalized();
        colision.collision_objects = {sphere_ind, id};

        return true;
        
    }
    else{
        //cout << " ----- " << endl;
        //cout << "No collision" << endl;
        colision.type = Collision::null;
        return false;   
    }

}

