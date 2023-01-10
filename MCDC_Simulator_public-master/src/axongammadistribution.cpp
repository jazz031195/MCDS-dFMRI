#include "axongammadistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;

AxonGammaDistribution::AxonGammaDistribution (double dyn_perc_,double volume_inc_perc_, unsigned num_ax, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l, float min_radius, bool active_state_)
{
    dyn_perc = dyn_perc_;
    volume_inc_perc = volume_inc_perc_;
    num_obstacles = num_ax;
    alpha = a;
    beta  = b;
    icvf = icvf_;
    min_limits = min_l;
    max_limits = max_l;
    axons.clear();
    this->min_radius = min_radius;
    active_state = active_state_;
    projections_x.clear();
    projections_y.clear();
    projections_z.clear();


    //for checking collision

}
void AxonGammaDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d &l)
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

void AxonGammaDistribution::displayGammaDistribution()
{
    const int nrolls=10000;  // number of experiments
    const int nstars=100;    // maximum number of stars to distribute
    string message;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha,beta);

    int p[11]={};

    for (int i=0; i<nrolls; ++i) {
        double number = distribution(generator);

        if (number<10) ++p[int(number)];
        else ++p[10];
    }

    for (int i=0; i<9; ++i) {
        message = std::to_string(i) + "-" + std::to_string(i+1) + ": " + std::string(p[i]*nstars/nrolls,'*');
        SimErrno::info(message,cout);
    }
    message = "9-10:" + std::string(p[9]*nstars/nrolls,'*') ;
    SimErrno::info(message,cout);
    message = ">10: " +  std::string(p[10]*nstars/nrolls,'*') + "\n" ;
    SimErrno::info(message,cout);
}


void AxonGammaDistribution::add_projection(Axon ax, int ax_index, double distance_to_be_inside, ostream& out){
    string message;
    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {

        for (int i=0; i<ax.spheres.size(); ++i){ 

            double position1, position2;
            int ax_id = ax_index;
            std::vector<projection_pt> axon_projection_on_axis;
            int sph_id = i;

            if (axis == 0){
                axon_projection_on_axis = projections_x;
            }
            else if (axis == 1){
                axon_projection_on_axis = projections_y;
            }
            else if (axis == 2){
                axon_projection_on_axis = projections_z;
            }

            // center + radius

            position1 = ax.spheres[i].center[axis] + ax.max_radius;
            //  find the first element smaller than the new position
            auto pos1 = std::find_if(axon_projection_on_axis.begin(), axon_projection_on_axis.end(), [position1](projection_pt s) {
                return s.position < position1 ;
            });
            // And then insert the new element at this position

            axon_projection_on_axis.insert(pos1, {position1, ax_id, sph_id});

            // center - radius

            position2 = ax.spheres[i].center[axis] - ax.max_radius;
            //  find the first element smaller than the new position
            auto pos2 = std::find_if(axon_projection_on_axis.begin(), axon_projection_on_axis.end(), [position2](projection_pt s) {
                return s.position < position2 ;
            });

            // And then insert the new element at this position
            axon_projection_on_axis.insert(pos2, {position2, ax_id, sph_id});
            
            if (axis == 0){
                projections_x = axon_projection_on_axis;
            }
            else if (axis == 1){
                projections_y = axon_projection_on_axis;
            }
            else if (axis == 2){
                projections_z = axon_projection_on_axis;
            }
        }
    }
    
} 
void AxonGammaDistribution::createGammaSubstrate(ostream& out)
{
    // generate the gamma distribution
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    uint repetition = 2;
    uint max_adjustments = 5;
    double best_icvf = 0;
    vector<Axon> best_axons;
    Eigen::Vector3d best_max_limits;
    min_limits = {0., 0., 0.};

    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0, 1);
    std::vector<double> radiis(num_obstacles, 0);

    int number_swelling_axons = int(num_obstacles * dyn_perc);
    //std::vector<int> swell_cyl_id(number_swelling_cylinders, 0);
    std::vector<bool> bool_swell_ax_id(num_obstacles, false);

    bool achieved = false;

    int tried = 0;

    string message;

    projections_x.clear();
    projections_y.clear();
    projections_z.clear();

    // create list of ids that correspond to dynamic axons based on percentage input
    for (unsigned i = 0; i < number_swelling_axons; ++i)
    {
        int random_id = rand() % num_obstacles;
        if (i > 0)
        {
            while (bool_swell_ax_id[random_id]) {
                random_id = rand() % num_obstacles;
            }
        }
        //swell_cyl_id[i] = random_id;
        bool_swell_ax_id[random_id] = true;
    }


    for (unsigned i = 0; i < num_obstacles; ++i)
    {
        

        if (tried > 10000)
        {
            message = " Radii distribution cannot be sampled [Min. radius Error]\n";
            SimErrno::error(message, cout);
            assert(0);
        }
        double jkr = distribution(generator);

        if (jkr < this->min_radius)
        {
            i--;
            tried++;
            continue;
        }
        tried = 0;

        radiis[i] = jkr * 1e-3; // WE CONVERT FROM UM TO MM HERE
    }


    // using a lambda function:
    std::sort(radiis.begin(), radiis.end(), [](const double a, double b) -> bool
              { return a > b; });

    double max_radius_ = radiis[0] * sqrt(1+volume_inc_perc);

    uint adjustments = 0;
    // We increease 1% the total area. (Is prefered to fit all the cylinders than achieve a perfect ICVF.)
    double adj_increase = icvf * 0.01;


    //while (!achieved)
    //{

        double target_icvf = this->icvf + adjustments * adj_increase;
        // computes max_limits
        computeMinimalSize(radiis, target_icvf, max_limits);


        for (uint t = 0; t < repetition; t++)
        {
            //vector<Axon> axons_to_add;

            axons.clear();
            for (unsigned i = 0; i < num_obstacles; i++)
            {

                //out << " obstacle :"<< i << endl;
                double t = udist(gen);
                double distance_to_border = radiis[i]*sqrt(1+volume_inc_perc) + barrier_tickness;
                double x = (t * (max_limits[0]-distance_to_border)) + (1 - t) * (min_limits[0]+ distance_to_border);
                t = udist(gen);
                double y = (t * (max_limits[1]-distance_to_border) + (1 - t) * (min_limits[1]+ distance_to_border));

                Vector3d Q = {x, y, min_limits[2] + distance_to_border};
                Vector3d D = {x, y, max_limits[2] - distance_to_border};

                Axon ax(radiis[i], Q, D, volume_inc_perc, active_state, bool_swell_ax_id[i], 1);

                std::vector<Dynamic_Sphere> spheres_to_add = GrowAxon(ax, max_radius_, i,  out);

                ax.set_spheres(spheres_to_add);


                //bool collision = isColliding(ax, max_radius + EPS_VAL, i, out);

                //if (!collision)
                if(spheres_to_add.size() != 0)
                {
                    axons.push_back(ax);
                    add_projection(ax, i, EPS_VAL, out);
                 }
                else{
                    i -= 1;
                    continue;
                }  

                int dummy;
                double icvf_current = computeICVF(axons, min_limits, max_limits, dummy);
                if (icvf_current > best_icvf)
                {
                    best_icvf = icvf_current;
                    best_axons.clear();
                    best_axons = axons;
                    best_max_limits = max_limits;
                }


            } // end for axons

            if (this->icvf - best_icvf < 0.0005)
            {
                achieved = true;
                break;
            }
        }
        axons.clear();
        adjustments++;
        cout << best_icvf << endl;
        if (adjustments > max_adjustments)
        {
            //break;
        }
    //}

    axons = best_axons;
    max_limits = best_max_limits;

    for (unsigned i = 0; i < axons.size(); i++){
        out << "Axon :" << i << endl;
        for (unsigned s = 0; s < axons[i].spheres.size(); s++){
            out << axons[i].spheres[s].center[0] << " " << axons[i].spheres[s].center[1] << " " << axons[i].spheres[s].center[2] << " " << axons[i].spheres[s].radius << endl;
        }
    }


    // TODO cambiar a INFO
    int perc_;
    double icvf_current = computeICVF(axons, min_limits, max_limits, perc_);

    message = "Percentage of axons selected: " + to_string(double(perc_) / radiis.size() * 100.0) + "%,\nICVF achieved: " + to_string(icvf_current * 100) + "  (" + to_string(int((icvf_current / icvf * 100))) + "% of the desired icvf)\n";
    SimErrno::info(message, cout);
    message = "number of axons :  " + to_string(axons.size()) +"\n";
    SimErrno::info(message, cout);

}

void AxonGammaDistribution::printSubstrate(ostream &out)
{
    out << 1e-3 << endl;
    out << volume_inc_perc << endl;
    out << dyn_perc << endl;
    out << icvf << endl;

    for (unsigned i = 0; i < axons.size(); i++)
    {
        out << axons[i].begin[0] * 1e3  << " " << axons[i].begin[1] * 1e3  << " " 
        << axons[i].begin[2] * 1e3  << " " << axons[i].end[0] * 1e3  << " " 
        << axons[i].end[1] * 1e3  <<" "<< axons[i].end[2] * 1e3  <<" " 
        << axons[i].min_radius* 1e3 << " " << axons[i].swell 
        << endl;

    }
}

std::vector<AxonGammaDistribution::projection_pt> AxonGammaDistribution::find_collisions(projection_pt proj_on_axis_min, projection_pt proj_on_axis_max,std::vector<projection_pt> projections_on_axis, ostream& out){
    
    std::vector<projection_pt> closest_spheres;
    closest_spheres.clear();
    string message;


    if (projections_on_axis.size()== 0){
        return closest_spheres; 
    } 

    // projection after which projections are smaller than min
    auto pos_min = std::find_if(projections_on_axis.begin(), projections_on_axis.end(), [proj_on_axis_min](projection_pt s) {
        return s.position <= proj_on_axis_min.position ;
    });

    // projection index after which projections are smaller than min
    unsigned index_min = std::distance(std::begin(projections_on_axis), pos_min);  

    // projection index after which projections are smaller than max
    auto pos_max = std::find_if(projections_on_axis.begin(), projections_on_axis.end(), [proj_on_axis_max](projection_pt s) {
        return s.position <= proj_on_axis_max.position ;
    });

    // projection index after which projections are smaller than max
    unsigned index_max = std::distance(std::begin(projections_on_axis), pos_max);

    if (index_min == index_max){
        return closest_spheres; 
    } 

    //out << " **---------------------find collision-------------------**" << endl;
    //out << "sphere : " << proj_on_axis_max.sphere_id << ", axon : " << proj_on_axis_max.axon_id << endl;
    //out << "proj_on_axis_max.position : " << proj_on_axis_max.position << endl;
    //out << "proj_on_axis_min.position : " << proj_on_axis_min.position << endl;
    //out << "index_min : " << index_min << endl;
    //out << "index_max :" << index_max << endl;

    //for (unsigned i = 0 ; i < projections_on_axis.size(); i++){
    //    out <<"pos : "<<projections_on_axis[i].position << " sphere :" << projections_on_axis[i].sphere_id << ", axon : "<< projections_on_axis[i].axon_id << endl;
    //    if (i == index_min){
    //        out << "--------- index_min ------------" << endl;
    //    } 
    //    if (i == index_max){
    //        out << "--------- index_max ------------" << endl;
    //    } 
    //} 


    for (unsigned i = index_max ; i < index_min; i++){
        
        projection_pt s{projections_on_axis[i].position, projections_on_axis[i].axon_id, projections_on_axis[i].sphere_id};  
        if (!search_for_sphere(closest_spheres,s) && s.axon_id != proj_on_axis_min.axon_id){ 
            closest_spheres.push_back(s);
        } 
        //out << "added : sphere_id : " <<s.sphere_id << ", axon_id : "<<s.axon_id <<  endl;

        //message = "close sphere : "+to_string(s.position)+"\n";
        //SimErrno::info(message,cout);
    } 

    return closest_spheres;
}

bool AxonGammaDistribution::search_for_sphere(std::vector<AxonGammaDistribution::projection_pt> spheres_, AxonGammaDistribution::projection_pt s){
    
    // search for s in spheres_
    if (spheres_.size() == 0){
        return false;
    } 
    for (unsigned i = 0; i < spheres_.size(); i++){
        if (spheres_[i].axon_id == s.axon_id && spheres_[i].sphere_id == s.sphere_id){
            return true;
            break;

        }  
    } 
    return false;
}  

bool AxonGammaDistribution::isSphereColliding(Dynamic_Sphere sph, double distance_to_be_inside, int axon_id, int sph_id, ostream& out){
    std::vector<std::vector<projection_pt>> coliding_spheres;
    std::vector<projection_pt> colisions_axis_spheres;
    coliding_spheres.clear();
    // on all axes
    for (unsigned x = 0; x < 3; x++){
        //out << "************" << x << "***********" << endl;
        colisions_axis_spheres.clear();

        projection_pt proj_on_axis_min = {sph.center[x]- distance_to_be_inside, axon_id, sph_id};
        // get max projection
        projection_pt proj_on_axis_max = {sph.center[x] + distance_to_be_inside, axon_id, sph_id};

        if (x== 0){

            colisions_axis_spheres = find_collisions(proj_on_axis_min, proj_on_axis_max, projections_x, out);

        }  
        else if (x == 1) {

            colisions_axis_spheres = find_collisions(proj_on_axis_min, proj_on_axis_max, projections_y, out);
        }  
        else{

            colisions_axis_spheres = find_collisions(proj_on_axis_min, proj_on_axis_max, projections_z, out);
            
        } 
        if (colisions_axis_spheres.size()== 0){
            return false;
            break;
        } 
        else{ 
            coliding_spheres.push_back(colisions_axis_spheres);

        } 
    }
    //out << "coliding_spheres.size : " << coliding_spheres.size() << endl;

    if (coliding_spheres.size() == 3){ 

        // for all coliding spheres in x 
        for(unsigned j = 0; j < coliding_spheres[0].size() ; j++){ 
            const projection_pt coliding_sphere = coliding_spheres[0][j];
                // if the same coliding spheres are also in y and z but are not from same axon 
            if (search_for_sphere(coliding_spheres[0], coliding_sphere) && search_for_sphere(coliding_spheres[1], coliding_sphere) && search_for_sphere(coliding_spheres[2], coliding_sphere)){
                
                //out << "coliding_sphere.position:" << coliding_sphere.position << ", axon_id : " << coliding_sphere.axon_id << ", sphere_id : "<< coliding_sphere.sphere_id << endl;;
                int sph_id = coliding_sphere.sphere_id;
                int ax_id = coliding_sphere.axon_id;  
                Dynamic_Sphere sphere_ = axons[ax_id].spheres[sph_id];

                //out << "sphere_ center :" << sphere_.center[0]  << ", center + radius : " << sphere_.center[0] + sphere_.radius  << ", center - radius : " << sphere_.center[0] - sphere_.radius << endl;
                //out << "axon radius :" << axons[ax_id].radius  << ", max_radius : " << axons[ax_id].max_radius   << ", swell : " << axons[ax_id].swell << endl;
                //out << "sphere radius" << sphere_.radius << endl;
                if (sph.isInside(sphere_.center, axons[ax_id].max_radius + 2*barrier_tickness)){ 
                    //out <<  "colision !" << endl;;
                    return true;
                    break;
                }  
            } 
        }
    }
    return false; 

} 
 
bool AxonGammaDistribution::isColliding(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out){
    string message;

    if (axons.size() == 0){
        return false;
    } 
    // for all spheres in axon
    for (unsigned i = 0; i < ax.spheres.size(); i++){
        out << "sphere :" << i << endl;
        bool colision = isSphereColliding(ax.spheres[i], distance_to_be_inside, axon_id, i, out);
        if (colision){
            return true;
            break;
        } 
    }
    return false;
}


/*
WARNING: The way we discard repeated cylinders is using radius. Repreated radius (like really the same)
are considered the same. This becasuse we don't track wich cylinders had to be replciated to mantain the voxel
symmetry
*/

double AxonGammaDistribution::computeICVF(std::vector<Axon> &axons, Vector3d &min_limits, Vector3d &max_limits, int &num_no_repeat)
{

    if (axons.size() == 0)
        return 0;

    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1])*(max_limits[2] - min_limits[2]);

    double AreaC = 0;

    // using a lambda function:
    std::sort(axons.begin(), axons.end(), [](const Axon a, Axon b) -> bool
              { return a.radius > b.radius; });

    double rad_holder = -1;
    num_no_repeat = 0;
    for (uint i = 0; i < axons.size(); i++)
    {

        if (fabs(rad_holder - axons[i].radius) < 1e-15)
        {
            continue;
        }
        else
        {
            rad_holder = axons[i].radius;
        }
        double rad = axons[i].radius;

        double ax_length = ((axons[i].spheres).size()-1)*rad/4;
        
        AreaC += M_PI * rad * rad* ax_length;
        num_no_repeat++;
    }
    return AreaC / AreaV;
}



std::tuple<double, double>  phi_gamma_to_target (Eigen::Vector3d prev_pos, Eigen::Vector3d new_pos, Eigen::Vector3d end,  ostream& out) {

    Eigen::Vector3d vector_to_target = {end[0] - new_pos[0], end[1] - new_pos[1], end[2] - new_pos[2]};
    Eigen::Vector2d vector_to_target_xy = {vector_to_target[0] ,vector_to_target[1]};
    vector_to_target_xy = vector_to_target_xy.normalized();
    Eigen::Vector2d vector_to_target_xz = {vector_to_target[0] ,vector_to_target[2]}; 
    vector_to_target_xz = vector_to_target_xz.normalized();

    //out << "vector_to_target_xy : ( " << vector_to_target_xy[0] << ", " << vector_to_target_xy[1] << " )" << endl;
    //out << "vector_to_target_xz : ( " << vector_to_target_xz[0] << ", " << vector_to_target_xz[1] << " )" << endl;

    Eigen::Vector3d straight_vector = {new_pos[0] - prev_pos[0], new_pos[1] - prev_pos[1], new_pos[2] - prev_pos[2]};
    Eigen::Vector2d straight_vector_xy = {straight_vector[0] ,straight_vector[1] };
    straight_vector_xy = straight_vector_xy.normalized();
    Eigen::Vector2d straight_vector_xz = {straight_vector[0] ,straight_vector[2] }; 
    straight_vector_xz = straight_vector_xz.normalized();

    //out << "straight_vector_xy : ( " << straight_vector_xy[0] << ", " << straight_vector_xy[1] << " )" << endl;
    //out << "straight_vector_xz : ( " << straight_vector_xz[0] << ", " << straight_vector_xz[1] << " )" << endl;


    double phi_to_target = atan2(vector_to_target_xy[1], vector_to_target_xy[0])-atan2(straight_vector_xy[1],straight_vector_xy[0]);
    double gamma_to_target = atan2(vector_to_target_xz[1],vector_to_target_xz[0])-atan2(straight_vector_xz[1],straight_vector_xz[0]);

    if (phi_to_target > M_PI)
    { phi_to_target -= 2 * M_PI; }
    else if (phi_to_target <= -M_PI)
    { phi_to_target += 2 * M_PI; }

    if (gamma_to_target > M_PI)
    { gamma_to_target -= 2 * M_PI; }
    else if (gamma_to_target <= -M_PI)
    { gamma_to_target += 2 * M_PI; }

    return std::make_tuple(phi_to_target, gamma_to_target);
}

bool AxonGammaDistribution::check_borders(Eigen::Vector3d pos, double distance_to_border) {
    Eigen::Vector3d new_min_limits = {min_limits[0] + distance_to_border, min_limits[1] + distance_to_border,min_limits[2] + distance_to_border};
    Eigen::Vector3d new_max_limits = {max_limits[0] - distance_to_border, max_limits[1] - distance_to_border,max_limits[2] - distance_to_border};

    if ((pos[0]-new_min_limits[0])<0 || (pos[1]-new_min_limits[1])<0 || (pos[2]-new_min_limits[2])<0){
        return true;
    }
    else if ((pos[0]-new_max_limits[0])>0 || (pos[1]-new_max_limits[1])>0 || (pos[2]-new_max_limits[2])>0) {
        return true;
    }
    return false;
    
}

std::vector<Dynamic_Sphere> AxonGammaDistribution::GrowAxon(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out){

    std::vector<Eigen::Vector3d> centers;
    double dist_ = ax.radius;
    Eigen::Vector3d prev_pos = ax.begin;
    centers.push_back(prev_pos);

    std::random_device rd;
    std::mt19937 gen(rd()); 

    int sphere_id = 0;

    double phi;
    double gamma;
    double phi_to_target, gamma_to_target;
    double delta_x;
    double delta_y;
    double delta_z;
    int tries ;
    bool achieved;
    Eigen::Vector3d pos_;
    Eigen::Vector3d new_pos = {prev_pos[0], prev_pos[1], prev_pos[2]+dist_};
    centers.push_back(new_pos);

    Dynamic_Sphere s0(prev_pos, ax.radius,ax.volume_inc_perc,ax.swell, axon_id, 1);
    Dynamic_Sphere s1(new_pos, ax.radius,ax.volume_inc_perc,ax.swell, axon_id, 1);
    std::vector<Dynamic_Sphere> spheres_to_add;

    int max_tries = 20;
    bool stop = false;

    out << "try_axon :" << axon_id << endl;


    if(isSphereColliding(s0, distance_to_be_inside, axon_id, sphere_id, out) || isSphereColliding(s1, distance_to_be_inside, axon_id, sphere_id, out)) {
        return spheres_to_add;
    }

    int fibre_collapse = 0;
    do{  

        //out << "pos : (" << new_pos[0] << ", " << new_pos[1] << ", " << new_pos[2] << ") " << endl;
        //out << "end: (" << ax.end[0] << ", " << ax.end[1] << ", " << ax.end[2] << ") " << endl;
        tries =0;
        tie(phi_to_target, gamma_to_target) = phi_gamma_to_target (prev_pos, new_pos, ax.end, out);
        out << "phi to target : "<< phi_to_target/M_PI << "*pi, gamma to target :"<< gamma_to_target/M_PI <<"*pi" << endl;
        achieved = false;

        while(!achieved && tries < max_tries){
            std::normal_distribution<float> phi_dist (phi_to_target/M_PI, 0.1); 
            phi = phi_dist(gen)*M_PI;
            std::normal_distribution<float> gamma_dist (gamma_to_target/M_PI, 0.1); 
            gamma = gamma_dist(gen)*M_PI;
            if (gamma > M_PI/4){
                gamma = M_PI/4;
            } 
            else if (gamma < - M_PI/4) {
                gamma = -M_PI/4;
            }

            out << "phi : " << phi/M_PI << "*pi , gamma : " << gamma/M_PI <<"*pi"<< endl;
            delta_x = dist_*cos(phi)*sin(gamma);
            delta_y = dist_*sin(phi)*sin(gamma);
            delta_z = dist_*cos(gamma);

            //out << "delta_y : "<< delta_y << ", delta_x :"<< delta_x <<", delta_z :" << delta_z<< endl;
        
            pos_ = new_pos;
            
            pos_[0] +=  delta_x;
            pos_[1] +=  delta_y;
            pos_[2] +=  delta_z;

            // check distance with border

            if (!check_borders(pos_, ax.max_radius + barrier_tickness)) {

                // check collision with other spheres
                Dynamic_Sphere s(pos_, ax.min_radius,ax.volume_inc_perc,ax.swell, axon_id, 1);

                if (!isSphereColliding(s, distance_to_be_inside, axon_id, sphere_id, out)){  
                    sphere_id += 4;
                    prev_pos = new_pos;
                    new_pos = pos_;
                    centers.push_back(new_pos);
                    achieved = true;
                } 
                else {
                    //out << "collision" << endl;
                    tries += 1;
                }
            }
            else {
                tries += 1;
            }
        }
        // fibre collapse
        if(tries >= max_tries && !achieved){
            if (centers.size() > 3 && fibre_collapse < 5){
                centers.pop_back();
                out << "Max tries reached ! centers.size() : " << centers.size() << endl;
                new_pos = centers[-1];
                prev_pos = centers[-2];
                fibre_collapse += 1;
            }
            else{
                return spheres_to_add;
            }
        }

    }while (new_pos[2] < ax.end[2]-dist_ );



    for (unsigned i=0; i< centers.size(); ++i){
        Dynamic_Sphere sphere(centers[i], ax.min_radius,ax.volume_inc_perc,  ax.swell, ax.id, 1);  
        spheres_to_add.push_back(sphere);
        if (i != 0){
            // squeleton has centers every r, so 3 spheres must be created in between

            Eigen::Vector3d center_between0 = {centers[i][0]/4+3*centers[i-1][0]/4, centers[i][1]/4+3*centers[i-1][1]/4, centers[i][2]/4+3*centers[i-1][2]/4};
            Dynamic_Sphere sphere0 (center_between0, ax.min_radius,ax.volume_inc_perc,  ax.swell, ax.id, 1);  
            spheres_to_add.push_back(sphere0);

            Eigen::Vector3d center_between1 = {(centers[i][0]+centers[i-1][0])/2, (centers[i][1]+centers[i-1][1])/2, (centers[i][2]+centers[i-1][2])/2};
            Dynamic_Sphere sphere1(center_between1, ax.min_radius,ax.volume_inc_perc,  ax.swell, ax.id, 1);  
            spheres_to_add.push_back(sphere1);

            Eigen::Vector3d center_between2 = {3*centers[i][0]/4+centers[i-1][0]/4, 3*centers[i][1]/4+centers[i-1][1]/4, 3*centers[i][2]/4+centers[i-1][2]/4};
            Dynamic_Sphere sphere2(center_between2, ax.min_radius,ax.volume_inc_perc,  ax.swell, ax.id, 1);  
            spheres_to_add.push_back(sphere2);
        //out << sphere.center[0] << " " << sphere.center[1] << " " << sphere.center[2] << " " << sphere.radius << endl;
        }
    }

    return spheres_to_add;

}   