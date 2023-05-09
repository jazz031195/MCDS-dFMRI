#include "axongammadistribution.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"
#include <chrono>



using namespace std;
using namespace Eigen;
using namespace std::chrono;

Growth::Growth (Axon* axon_to_grow_, std::vector<Axon> env_axons_, Eigen::Vector3d small_voxel_size_, bool tortuous_){
           
    twin_axons.clear();
    env_axons = env_axons_;
    axon_to_grow = axon_to_grow_;
    axon_to_grow->id  = env_axons_.size();
    small_voxel_size = small_voxel_size_;
    tortuous = tortuous_;
    success = false;
    finished = false;
}

Growth::~Growth() { delete axon_to_grow; }


void Growth::add_next_sphere(Dynamic_Sphere added_sphere, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii){
    
    centers.push_back(added_sphere.center);
    sph_radii.push_back(added_sphere.radius);
        
}

void Growth::shrink_sphere_rad(double& rad, double axon_rad, double& shrink_perc){
    rad = axon_rad *(1-shrink_perc);

}

void Growth::find_shrinking_dichotomy(double& rad, double axon_rad, double min_shrink, double max_shrink, double& half_shrink){

    half_shrink = (max_shrink-min_shrink)/2+min_shrink;
    rad = axon_rad *(1-half_shrink);
    if (rad < 0.05*1e-3) {
        rad = 0.05*1e-3;
    }

}
bool Growth::check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d& twin_delta_pos) {
    
    
    Eigen::Vector2d new_min_limits = {distance_to_border, distance_to_border};
    Eigen::Vector2d new_max_limits = {small_voxel_size[0] - distance_to_border, small_voxel_size[1] - distance_to_border};

    twin_delta_pos = {0.0,0.0};

    if ((pos[0]-new_min_limits[0])<0) {
        // min plane of x
        twin_delta_pos[0] = small_voxel_size[0];
    }
    if ((pos[1]-new_min_limits[1])<0 ){
        // min plane of y
        twin_delta_pos[1] = small_voxel_size[1];
    } 
    if ((pos[0]-new_max_limits[0])>0) {
        // max plane of x
        twin_delta_pos[0] = -small_voxel_size[0];
    }
    if ((pos[1]-new_max_limits[1])>0 ){
        // max plane of y
        twin_delta_pos[1] = -small_voxel_size[1];
    } 
    if (twin_delta_pos !=  Eigen::Vector2d({0.0,0.0})){
        return true;
    }   
    else {
        return false;
    } 
    
}

std::tuple<double, double>  phi_theta_to_target (Eigen::Vector3d new_pos, Eigen::Vector3d end,  ostream& out) {
    /*
    Find phi and theta angles in spherical coordinates
    see : http://douglashalse.com/index.php/2019/09/25/spherical-coordinates/

    */

    Eigen::Vector3d vector_to_target;
    if (end[2]< new_pos[2]){ 
        vector_to_target= {end[0] - new_pos[0], end[1] - new_pos[1], end[2]-0.01 - new_pos[2]};
    }
    else{
        vector_to_target= {end[0] - new_pos[0], end[1] - new_pos[1], end[2]+0.01 - new_pos[2]};
    }  
    vector_to_target = vector_to_target.normalized();

    double phi_to_target;
    double theta_to_target;

    // theta is angle between (1,0) and (x,y)
    // varies between 0 and 2pi
    theta_to_target = atan2(vector_to_target[1],vector_to_target[0]);
    if (theta_to_target < 0){
        theta_to_target += 2*M_PI;
    } 
    
    // phi is angle between (0,0,1) and (x,y,z)
    // varies between 0 and pi
    if (vector_to_target[2] == 0){
        phi_to_target = M_PI/2;
    } 
    else if (vector_to_target == Eigen::Vector3d({0,0,-1})){
        phi_to_target = M_PI;
    } 
    else if (vector_to_target == Eigen::Vector3d({0,0,1})){
        phi_to_target = 0;
    } 
    else{ 
        // varies between -pi/2 and pi/2 
        phi_to_target = atan((sqrt(vector_to_target[0]*vector_to_target[0]+vector_to_target[1]*vector_to_target[1]))/vector_to_target[2]);
        if (phi_to_target < 0){
            phi_to_target += M_PI;
        } 
    } 

    return std::make_tuple(phi_to_target, theta_to_target);
}

bool Growth::isSphereColliding(Dynamic_Sphere sph){
    
    Vector3d position = sph.center;

    std::vector<int> col_sphere_ids;

    //std::cout << "env_axons.size():" << env_axons.size() << endl;

    for (unsigned i = 0; i < env_axons.size() ; i++){
        double distance_to_be_inside = sph.radius;
        bool isinside = env_axons[i].isPosInsideAxon(position, distance_to_be_inside, col_sphere_ids);
        
        if (isinside){
            return true;
            break;
        }
    }
    return false;
}

void Growth::createTwinSpheres(std::vector<Dynamic_Sphere>& twin_spheres, Dynamic_Sphere s, Eigen::Vector2d twin_delta_pos){

    Dynamic_Sphere twin_sphere = s;
    twin_sphere.center = {twin_sphere.center[0]+twin_delta_pos[0], twin_sphere.center[1]+twin_delta_pos[1], twin_sphere.center[2]};
    twin_spheres.push_back(twin_sphere);

    if (twin_delta_pos[0] != 0 && twin_delta_pos[1]!= 0){
        Dynamic_Sphere twin_sphere1 = s;
        twin_sphere1.center = {twin_sphere1.center[0], twin_sphere1.center[1]+twin_delta_pos[1], twin_sphere1.center[2]};
        twin_spheres.push_back(twin_sphere1);

        Dynamic_Sphere twin_sphere2 = s;
        twin_sphere2.center = {twin_sphere2.center[0]+twin_delta_pos[0], twin_sphere2.center[1], twin_sphere2.center[2]};
        twin_spheres.push_back(twin_sphere2);
    }
}

bool Growth::find_next_center(Axon* ax, Eigen::Vector3d destination, Dynamic_Sphere& s, vector<Eigen::Vector3d> centers, double dist_, double& rad, ostream& out, bool& collides_with_border, Eigen::Vector2d& twin_delta_pos,  std::vector<Eigen::Vector2d> phi_theta_colliding, int max_tries = 300){

    // phi and theta angles in spherical coordinates
    // to update position with
    double phi, theta;
    // phi and theta angles in spherical coordinates
    // that links the current position to the target position
    double phi_to_target, theta_to_target;
    // phi and theta angles in spherical coordinates
    // that links previous position with current one
    double prev_phi, prev_theta;
    // number of tries 
    int tries = 0;
    // random
    std::random_device rd;
    std::mt19937 gen(rd()); 

    // differences in position from intial to final position
    double delta_x;
    double delta_y;
    double delta_z;

    // position of new sphere to add
    Eigen::Vector3d new_pos;
    // current position (position of last sphere)
    Eigen::Vector3d curr_pos = {centers[centers.size() -1][0],centers[centers.size() -1][1],centers[centers.size() -1][2]};
    // current position (position of 4 spheres before last )
    Eigen::Vector3d prev_pos;
    if (centers.size()>2){
        prev_pos = {centers[centers.size()-2][0],centers[centers.size() -2][1],centers[centers.size() -2][2]};
    }
    else {
        prev_pos = {centers[centers.size()-1][0],centers[centers.size() -1][1],centers[centers.size() -1][2]-dist_};
    }

    while(tries < max_tries){

        // if collides with borders
        collides_with_border = false;

        tie(phi_to_target, theta_to_target) = phi_theta_to_target (curr_pos, destination, out);

        // if not tortuous, phi and theta are those that lead to target
        if (!tortuous){
            phi = phi_to_target;
            theta = theta_to_target;
        }
        // if tortuous, phi and theta comme from distribution
        else{
            std::normal_distribution<float> phi_dist (phi_to_target/M_PI, 0.08); 
            std::normal_distribution<float> theta_dist (theta_to_target/M_PI, 0.08); 
            phi = phi_dist(gen)*M_PI;

            
            // if phi has already been tried and didn't work in a combination
            // find all indexes at which the chosen phi has already been tried and hasn't worked
            std::vector<size_t> exclude_ind;
            double delta = 0.05;
            auto it = std::find_if(std::begin(phi_theta_colliding), std::end(phi_theta_colliding), [delta, phi](Eigen::Vector2d i){return (i[0] > phi/M_PI-delta && i[0] < phi/M_PI+delta);});
            while (it != std::end(phi_theta_colliding)) {
                exclude_ind.emplace_back(std::distance(std::begin(phi_theta_colliding), it));
                it = std::find_if(std::next(it), std::end(phi_theta_colliding), [delta, phi](Eigen::Vector2d i){return (i[0]  > phi/M_PI-delta && i[0] < phi/M_PI+delta);});
            }

            if (exclude_ind.size()>0){
                vector<double> theta_to_exclude;
                // find all the thetas to exclude because they were paired with phi
                for (unsigned i = 0; i < exclude_ind.size() ; i++){
                    theta_to_exclude.push_back(phi_theta_colliding[exclude_ind[i]][1]);
                } 
                // generate theta until it is not longer a theta that was paired with a phi that didn't work
                theta = theta_dist(gen)*M_PI;
                while (std::all_of(theta_to_exclude.cbegin(), theta_to_exclude.cend(), [theta, delta](double i){ return (i > theta/M_PI-delta && i< theta/M_PI+delta);})){
                    theta = theta_dist(gen)*M_PI;
                } 
            } 
            else{
                theta = theta_dist(gen)*M_PI;
            }  

        }

        // spherical coordinates to cartesian
        delta_x = dist_*cos(theta)*sin(phi);
        delta_y = dist_*sin(theta)*sin(phi);
        delta_z = dist_*cos(phi);
        

        // update new_pos
        new_pos = curr_pos;
        new_pos[0] +=  delta_x;
        new_pos[1] +=  delta_y;
        new_pos[2] +=  delta_z;

        // if angle difference is too much
        
        Eigen::Vector3d v1 = (curr_pos - prev_pos).normalized();
        Eigen::Vector3d v2 = (new_pos - curr_pos).normalized();
        

        Dynamic_Sphere sphere_ (new_pos, rad,ax->volume_inc_perc,ax->swell, ax->id, -1);
        
        // if sphere doesn't collide with environment
        if (!isSphereColliding(sphere_)){

            // if sphere collides with border planes
            if (check_borders(sphere_.center, sphere_.radius, twin_delta_pos)){
                std::vector<Dynamic_Sphere> twin_spheres;
                collides_with_border = true;
                createTwinSpheres(twin_spheres, sphere_, twin_delta_pos);
                int nbr_non_colliding_twins = 0;
                for (unsigned j=0; j< twin_spheres.size(); ++j){
                    if(!isSphereColliding(twin_spheres[j])){
                        nbr_non_colliding_twins +=1;

                    }
                }
                // if all twins don't collide with environment
                if (nbr_non_colliding_twins == twin_spheres.size()){
                    s = sphere_; 
                    return true;
                    break;
                }
                else {
                    tries += 1;
                    phi_theta_colliding.push_back({phi, theta});
                }
            }
            // if sphere doesn't collide with border planes
            else {
                collides_with_border = false;
                s =sphere_; 
                return true;
                break;

            }
        }
        else {
            tries += 1;
            phi_theta_colliding.push_back({phi, theta});
        }
    }
    
    return false;
}

bool Growth::fiber_collapse(std::vector<Eigen::Vector3d>& centers, int fibre_collapsed_nbr, ostream& out){
    
    int fiber_collapse_threshold = 5;
    // we delete 4 spheres at a time
    int nbr_discarded_centers = 1;
    if ((centers.size()> nbr_discarded_centers+1) && (fibre_collapsed_nbr <fiber_collapse_threshold)){
        for (unsigned j=0; j< nbr_discarded_centers; ++j){
            centers.pop_back();
        }
        return true;
    }
    else{
        return false;
    }
}



void Growth::fill_with_spheres(Axon* ax, std::vector<Dynamic_Sphere>& spheres_to_add, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii){

    for (unsigned i = 1; i < centers.size(); i++){

        Dynamic_Sphere sphere (centers[i-1], sph_radii[i-1] ,ax->volume_inc_perc, ax->swell, ax->id, 1, (i-1)*4);  
        spheres_to_add.push_back(sphere);
        // add spheres in between 
        for (unsigned j = 0; j < 3; j++){
            // center
            Eigen::Vector3d center_between = {(3-j)*centers[i-1][0]/4+(j+1)*centers[i][0]/4, (3-j)*centers[i-1][1]/4+(j+1)*centers[i][1]/4, (3-j)*centers[i-1][2]/4+(j+1)*centers[i][2]/4};
            // radius
            double radius_in_between = (3-j)*sph_radii[i-1]/4+(j+1)*sph_radii[i]/4;
            // sphere
            Dynamic_Sphere sphere_in_between (center_between, radius_in_between,ax->volume_inc_perc, ax->swell, ax->id, 1, (i-1)*4+j+1);  
            double shrink_perc = 0;
            Dynamic_Sphere sphere_ = sphere_in_between;
            // check if this sphere collides
            // at some point it should always find a radius at which it doesn't collide with environment
            // because it is between two overlapping spheres that do not collide with environment
            while (isSphereColliding(sphere_)){
                double rad;
                if (shrink_perc >1){
                    cout << "Error : cannot shrink spheres in between so that they don't collide with environment !" << endl;
                } 
                shrink_sphere_rad(rad, radius_in_between, shrink_perc);
                sphere_ = Dynamic_Sphere (center_between, rad,ax->volume_inc_perc, ax->swell, ax->id, 1, sphere_.id);  
                shrink_perc += 0.01;
            }
            sphere_in_between = sphere_;
            spheres_to_add.push_back(sphere_in_between);
            /*
            
            bool last_sphere_add =false;
            if (!last_sphere_add){  
                if (sphere_in_between.center[2]<= small_voxel_size[2]){
                    spheres_to_add.push_back(sphere_in_between);
    
                    if (i == centers.size()-1 && j== 2 ){
                        Dynamic_Sphere sphere (centers[i], sph_radii[i] ,ax.volume_inc_perc, ax.swell, ax.id, 1, ax.active_state, (i-1)*4+j+2);  
                        spheres_to_add.push_back(sphere);

                    } 
                }
                else{

                        // vector from centers[i-1] to centers[i] 
                        Eigen::Vector3d vector = (centers[i] - centers[i-1]).normalized();
                        // centers[i-1][2] + lambda*vector[2] = max_limits[2]
                        double lambda_ = (small_voxel_size[2]-centers[i-1][2])/vector[2];
                        Eigen::Vector3d last_center = centers[i-1] +lambda_*vector;
                        double weight_i_1 = abs((last_center-centers[i]).norm())/abs((centers[i-1]-centers[i]).norm());
                        double weight_i = abs((last_center-centers[i-1]).norm())/abs((centers[i-1]-centers[i]).norm());
                        double last_radius = weight_i*sph_radii[i]+weight_i_1*sph_radii[i-1];
                        Dynamic_Sphere last_sphere(last_center, last_radius,ax.volume_inc_perc,ax.swell, ax.id, 1, ax.active_state, (i-1)*4+j+1);
                        last_sphere_add = true;
                        if (!isSphereColliding(last_sphere)){
                            spheres_to_add.push_back(last_sphere);
                        }
                }
            } 
            */
        }
        if (i == centers.size()-1){
            Dynamic_Sphere sphere (centers[i], sph_radii[i] ,ax->volume_inc_perc, ax->swell, ax->id, 1, i*4);  
            spheres_to_add.push_back(sphere);

        }  
    } 
} 



std::vector<Dynamic_Sphere> Growth::GrowAxon(Axon* ax, Eigen::Vector3d destination, std::vector<Eigen::Vector2d>& all_twin_delta_pos, ostream& out){
    
    // list of all positions of the spheres add to axon in the growth
    std::vector<Eigen::Vector3d> centers;
    // radii of all spheres add to the centers list
    std::vector<double> sph_radii;
    // list of spheres to add
    std::vector<Dynamic_Sphere> spheres_to_add;

    // maximum of shrinking percentage
    double max_shrinking = 0.5;
    double min_radius_shrinking = ax->radius*(1-max_shrinking);
    // minimum radius of axon
    if (min_radius_shrinking < 0.05*1e-3) {
        min_radius_shrinking = 0.05*1e-3;
    }
    // managed to add 4 spheres
    bool spheres_added = false;
    // fiber collapsing is possible
    bool can_fiber_collapse = true;
    int fibre_collapsed_nbr = 0;
    // phi and thetas to exclude in tries
    std::vector<std::vector<Eigen::Vector2d>> all_phi_theta_colliding;


    // first sphere to add
    Dynamic_Sphere s1(ax->begin, ax->min_radius,ax->volume_inc_perc,ax->swell, ax->id, 1);

    
    if(Growth::isSphereColliding(s1)) {
        //std::cout << "first sphere collides "<< endl;
        return spheres_to_add;
    }
    else{
        // if the first sphere touches boundaries
        Eigen::Vector2d twin_delta_pos_;
        if(check_borders(s1.center, s1.radius, twin_delta_pos_)){
            std::vector<Dynamic_Sphere> first_twin_spheres;
            createTwinSpheres(first_twin_spheres, s1, twin_delta_pos_);
            int non_colliders = 0;
            // check if all twins don't collide with environment 
            for (unsigned i = 1; i < first_twin_spheres.size(); i++){
                if (!isSphereColliding(first_twin_spheres[i])){
                    non_colliders += 1;
                }
            }
            if (non_colliders == first_twin_spheres.size()) {
                add_next_sphere(s1, centers, sph_radii);
            }
            else{
                //std::cout << "twin of first sphere collides "<< endl;
                return spheres_to_add;
            }

        }
        else{
            add_next_sphere(s1, centers, sph_radii);
        }
    }
    bool stop_condition = false;

    do{
        // difference in position in case sphere goes out of boundaries
        Eigen::Vector2d twin_delta_pos;

        bool collides_with_border;

        double shrink_perc = 0.01;
        // radius of sphere to add
        double rad = ax->min_radius;
        // initial distance between spheres is the radius
        double dist_ = rad;

        Dynamic_Sphere sphere_to_add;
        std::vector<Eigen::Vector2d> phi_theta_colliding;
        // find 4 spheres to add that don't collide with environment at dist_ specified
        spheres_added = find_next_center(ax, destination, sphere_to_add, centers, dist_, rad, out, collides_with_border, twin_delta_pos, phi_theta_colliding);
        
        if (spheres_added){
            add_next_sphere(sphere_to_add, centers, sph_radii);
            // phi and theta that didn't work are saved 
            all_phi_theta_colliding.push_back(phi_theta_colliding);
            if(collides_with_border){

                // if all_twin_delta_pos is empty or twin_delta_pos not already in all_twin_delta_pos
                if (all_twin_delta_pos.size() == 0 || std::find(all_twin_delta_pos.begin(), all_twin_delta_pos.end(),twin_delta_pos) ==all_twin_delta_pos.end()) {
                    all_twin_delta_pos.push_back(twin_delta_pos);
                }
            }

            //out << "managed to add spheres :" << spheres_added<< endl;
            //out << "centers size :" << centers.size() << endl;
            //out << "all_phi_theta_colliding size :" << all_phi_theta_colliding.size() << endl;
        }

        else{

            // fiber collapse 
            // check and can do fiber collapsing
            can_fiber_collapse = fiber_collapse(centers, fibre_collapsed_nbr, out);
            if (can_fiber_collapse && all_phi_theta_colliding.size() > 0){
                std::vector<Eigen::Vector2d> phi_theta_colliding_collapse;
                //out << "fiber collapse"<< endl;
                // delete all necessary phi and thetas from previous steps
                all_phi_theta_colliding.pop_back();
                // find next center while not retaking phis and thetas that didn't work
                phi_theta_colliding_collapse = all_phi_theta_colliding[all_phi_theta_colliding.size()-1]; 
                // add next center
                spheres_added = find_next_center(ax, destination, sphere_to_add, centers, dist_, rad, out, collides_with_border, twin_delta_pos, phi_theta_colliding_collapse);
                if (spheres_added){
                    add_next_sphere(sphere_to_add, centers, sph_radii);
                    // add thetas and phis that didn't work in list
                    all_phi_theta_colliding.push_back(phi_theta_colliding_collapse);
                    // border
                    if(collides_with_border){
                        if (all_twin_delta_pos.size() == 0 || std::find(all_twin_delta_pos.begin(), all_twin_delta_pos.end(),twin_delta_pos) ==all_twin_delta_pos.end()) {
                            all_twin_delta_pos.push_back(twin_delta_pos);
                        }
                    }
                }
                fibre_collapsed_nbr += 1;
                //out << "centers size :" << centers.size() << endl;

            }
            // if cannot fiber collapse -> shrinking 
            if (!spheres_added){ 
                std::vector<Eigen::Vector2d> phi_theta_colliding_shrink;
                Dynamic_Sphere sphere_to_add_;
                int max_tries = 500;
                phi_theta_colliding_shrink.clear();
                bool can_add_with_max_shrinking = find_next_center(ax, destination, sphere_to_add_, centers, dist_, min_radius_shrinking, out, collides_with_border, twin_delta_pos, phi_theta_colliding_shrink, max_tries);
                // if a sphere can be added after a maximum shrinking process
                if (can_add_with_max_shrinking){
                    double lower_boundary = shrink_perc;
                    double upper_boundary = max_shrinking;
                    double middle_boundary;
                    bool dichotomy_check;
                    bool achieved = false;
                    Dynamic_Sphere last_sphere_to_add_= sphere_to_add_;
                    // while spheres can't be added and the radius doesn't reach the minimum required
                    while (!achieved){
                        // if boundaries are close enough
                        if (abs(upper_boundary-lower_boundary) < 0.01){
                            // shrink radius to the upper bound
                            add_next_sphere(last_sphere_to_add_, centers, sph_radii);
                            // add thetas and phis thta didn't work to remember
                            all_phi_theta_colliding.push_back(phi_theta_colliding_shrink);
                            if(collides_with_border){

                                // if all_twin_delta_pos is empty or twin_delta_pos not already in all_twin_delta_pos
                                if (all_twin_delta_pos.size() == 0 || std::find(all_twin_delta_pos.begin(), all_twin_delta_pos.end(),twin_delta_pos) ==all_twin_delta_pos.end()) {
                                    all_twin_delta_pos.push_back(twin_delta_pos);
                                }
                            }
                            //out << "radius shrink to :" << last_sphere_to_add_.radius << endl;
                            //out << "centers size :" << centers.size() << endl;
                            achieved = true;
                            }
                        else{   
                            // rad = shrinked radius
                            find_shrinking_dichotomy(rad, ax->radius, lower_boundary, upper_boundary, middle_boundary);
                            // distance between spheres is the maximum between the radius of the two
                            dist_ = max(rad, sph_radii[sph_radii.size()-1]);
                            // try shrinking at half 
                            phi_theta_colliding_shrink.clear();
                            dichotomy_check = find_next_center(ax, destination, sphere_to_add_, centers, dist_, rad, out, collides_with_border, twin_delta_pos, phi_theta_colliding_shrink);
                                
                            if (dichotomy_check){
                                upper_boundary = middle_boundary;
                                last_sphere_to_add_ = sphere_to_add_;
                            }
                            else{
                                lower_boundary = middle_boundary;
                            }
                        } 

                    }
                }
                else{
                    //out << "couldn't shrink more" << endl;
                    return spheres_to_add;
                    break;

                }  

            }
                
        }
    if (centers.size()>0) {
        if (destination[2]> s1.center[2] ){
            stop_condition= centers[centers.size()-1][2] > destination[2];
        } 
        else{
            stop_condition= centers[centers.size()-1][2] < destination[2];
        } 
    }
    else{
        stop_condition= false;
    }
    //out << "stop condition : " << stop_condition << endl;
    //out << "position : " << centers[centers.size()-1][2] << endl;

    //out << "small voxel size :" << small_voxel_size << endl;
    }while (!stop_condition);
    
    fill_with_spheres(ax, spheres_to_add,centers, sph_radii);

    return spheres_to_add;

}   

const void Growth::GrowInParallel(const Eigen::Vector3d destination){


    std::vector<Eigen::Vector2d> all_twin_delta_pos;


    std::vector<Dynamic_Sphere> spheres_to_add = GrowAxon(axon_to_grow, destination, all_twin_delta_pos, std::cout);
    //out <<"set spheres" << endl;
    if(spheres_to_add.size() != 0)
    {
        // sort by z comp of center
        std::sort(spheres_to_add.begin(),spheres_to_add.end(), [](const Dynamic_Sphere a, Dynamic_Sphere b) -> bool
              { return a.center[2] < b.center[2]; });
        axon_to_grow->set_spheres(spheres_to_add);
        success = true;
    }
    //out <<"create twins" << endl;
    // add twins if any
    for (unsigned k = 0; k < all_twin_delta_pos.size(); k++){
        Eigen::Vector2d twin_delta_pos = all_twin_delta_pos[k];
        createTwinAxons(axon_to_grow, twin_delta_pos);
    }
    finished = true;

}

void Growth::createTwinAxons(Axon* ax, Eigen::Vector2d twin_delta_pos){
    
    int id = env_axons.size()+twin_axons.size();
    Axon* twin_ax = new Axon (id, ax->min_radius, ax->begin, ax->end, ax->volume_inc_perc,  ax->swell, -1);
    std::vector<Dynamic_Sphere> twin_spheres;

    for (unsigned i = 0; i < ax->spheres.size() ; i++){
        Dynamic_Sphere twin_sphere = ax->spheres[i];
        twin_sphere.center = {twin_sphere.center[0]+twin_delta_pos[0], twin_sphere.center[1]+twin_delta_pos[1], twin_sphere.center[2]};
        twin_spheres.push_back(twin_sphere);
    }
    twin_ax->set_spheres(twin_spheres);
    twin_axons.push_back(*twin_ax);

    delete twin_ax;

    // if axon overlaps with 2 different planes of boundary
    if (twin_delta_pos[0] != 0 && twin_delta_pos[1] != 0){
        Axon* twin_ax1 = new Axon (id+1, ax->min_radius, ax->begin, ax->end, ax->volume_inc_perc, ax->swell, -1);
        twin_spheres.clear();

        for (unsigned i = 0; i < ax->spheres.size() ; i++){
            Dynamic_Sphere twin_sphere = ax->spheres[i];
            twin_sphere.center = {twin_sphere.center[0], twin_sphere.center[1]+twin_delta_pos[1], twin_sphere.center[2]};
            twin_spheres.push_back(twin_sphere);
        }
        twin_ax1->set_spheres(twin_spheres);
        twin_axons.push_back(*twin_ax1);
        delete twin_ax1;

        Axon* twin_ax2 = new Axon (id+2, ax->min_radius, ax->begin, ax->end, ax->volume_inc_perc,  ax->swell, -1);
        twin_spheres.clear();
        for (unsigned i = 0; i < ax->spheres.size() ; i++){
            Dynamic_Sphere twin_sphere = ax->spheres[i];
            twin_sphere.center = {twin_sphere.center[0]+twin_delta_pos[0], twin_sphere.center[1], twin_sphere.center[2]};
            twin_spheres.push_back(twin_sphere);
        }
        twin_ax2->set_spheres(twin_spheres);
        twin_axons.push_back(*twin_ax2);
        delete twin_ax2;
    }
}