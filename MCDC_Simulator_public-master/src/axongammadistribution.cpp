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
    icvf_current = 0;

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




void AxonGammaDistribution::createGammaSubstrate(ostream& out)
{
    // generate the gamma distribution
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    uint repetition = 1;
    uint max_adjustments = 2;
    double best_icvf = 0;
    vector<Axon> best_axons;
    Projections best_projections;
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
        

        if (tried > 1000)
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

            unsigned stuck = 0;

            for (unsigned i = 0; i < num_obstacles; i++)
            {
                if(stuck < 3000){
                    double t = udist(gen);
                    double distance_to_border = radiis[i]*sqrt(1+volume_inc_perc) + barrier_tickness;
                    double x = (t * (max_limits[0]-distance_to_border)) + (1 - t) * (min_limits[0]+ distance_to_border);
                    t = udist(gen);
                    double y = (t * (max_limits[1]-distance_to_border) + (1 - t) * (min_limits[1]+ distance_to_border));

                    //out << " obstacle :"<< i << endl;

                    Vector3d Q = {x, y, min_limits[2]};
                    Vector3d D = {x, y, max_limits[2]};


                    Axon ax(radiis[i], Q, D, volume_inc_perc, active_state, bool_swell_ax_id[i], 1);

                    std::vector<Dynamic_Sphere> spheres_to_add = GrowAxon(ax, max_radius_, i,  out);

                    if(spheres_to_add.size() != 0)
                    {
                        ax.set_spheres(spheres_to_add, i);
                        axons.push_back(ax);
                        stuck = 0;
                    }
                    else{
                        i -= 1; 
                        stuck += 1;
                        continue;
                    }  
                }

                int dummy;
                icvf_current = computeICVF(axons, min_limits, max_limits, dummy);
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
                //break;
            }
        //}
        axons.clear();
        adjustments++;
        cout << best_icvf << endl;
        if (adjustments > max_adjustments)
        {
            break;
        }
    }

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
    icvf_current = computeICVF(axons, min_limits, max_limits, perc_);

    out <<"icvf:"<< icvf_current << "voxel size: "<< max_limits[0] << endl;

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
    out << icvf_current << endl;
    out << min_limits[2] << endl;
    out << max_limits[2] << endl;
    
    for (unsigned i = 0; i < axons.size(); i++)
    {
        for (unsigned s = 0; s < axons[i].spheres.size(); s++){
            out << axons[i].spheres[s].center[0]*1e3 << " " << axons[i].spheres[s].center[1]*1e3 << " " << axons[i].spheres[s].center[2]*1e3 << " " << axons[i].spheres[s].radius*1e3 << " "<<axons[i].spheres[s].swell << endl;
        }
        out << "Axon: " << i << endl;
    }
}

double AxonGammaDistribution::computeICVF(std::vector<Axon> &axons, Vector3d &min_limits, Vector3d &max_limits, int &num_no_repeat)
{

    if (axons.size() == 0)
        return 0;

    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * (max_limits[2] - min_limits[2]);

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

        double ax_length = ((axons[i].spheres).size()-2)*rad/4;
        
        AreaC += M_PI * rad * rad* ax_length;
        num_no_repeat++;
    }
    return AreaC / AreaV;
}



std::tuple<double, double, double, double>  phi_gamma_to_target (Eigen::Vector3d prev_pos, Eigen::Vector3d new_pos, Eigen::Vector3d end,  ostream& out) {

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


    double phi_to_target = atan2(vector_to_target_xy[1], vector_to_target_xy[0]);
    double gamma_to_target = atan2(vector_to_target_xz[0],vector_to_target_xz[1]);

    // caculate angle wrt x axis
    double phi_straight = atan2(straight_vector_xy[1],straight_vector_xy[0]);
    // calculate angle wrt z axis
    double gamma_straight = atan2(straight_vector_xz[0],straight_vector_xz[1]);

    return std::make_tuple(phi_to_target, gamma_to_target, phi_straight, gamma_straight);
}

bool AxonGammaDistribution::check_borders(Eigen::Vector3d pos, double distance_to_border) {
    Eigen::Vector3d new_min_limits = {min_limits[0] + distance_to_border, min_limits[1] + distance_to_border,min_limits[2] + distance_to_border};
    Eigen::Vector3d new_max_limits = {max_limits[0] - distance_to_border, max_limits[1] - distance_to_border,max_limits[2] - distance_to_border};

    if ((pos[0]-new_min_limits[0])<0 || (pos[1]-new_min_limits[1])<0 ){
        return true;
    }
    else if ((pos[0]-new_max_limits[0])>0 || (pos[1]-new_max_limits[1])>0 ) {
        return true;
    }
    return false;
    
}

bool AxonGammaDistribution::isSphereColliding(Dynamic_Sphere sph, ostream& out){
    
    Vector3d position = sph.center;
    double distance_to_be_inside = 2* barrier_tickness + sph.max_radius;
    int ax_id;

    for (unsigned i = 0; i < axons.size() ; i++){
        bool isinside = axons[i].isPosInsideAxon(position, distance_to_be_inside, true);
        if (isinside){
            return true;
            break;
        }

    }
    return false;
}

bool AxonGammaDistribution::find_next_center(Axon ax, std::vector<double>& distances, double dist_,std::vector<double>& sph_radii, double& rad, Eigen::Vector3d& new_pos, Eigen::Vector3d& prev_pos, std::vector<Eigen::Vector3d>& centers, int axon_id, ostream& out){

    double phi;
    double gamma;
    double phi_to_target, gamma_to_target, gamma_straight, phi_straight;
    bool achieved = false;
    int tries = 0;
    int max_tries = 10000;
    int sphere_id = 0;
    Eigen::Vector3d pos_;

    std::random_device rd;
    std::mt19937 gen(rd()); 

    double delta_x;
    double delta_y;
    double delta_z;

    while(!achieved && tries < max_tries){
        std::normal_distribution<float> phi_dist (phi_to_target/M_PI, 0.08); 
        phi = phi_dist(gen)*M_PI;
        //phi = phi_to_target;
        std::normal_distribution<float> gamma_dist (gamma_to_target/M_PI, 0.08); 
        gamma = gamma_dist(gen)*M_PI;
        //gamma = gamma_to_target;
        if (gamma > M_PI/4+ gamma_straight){
            gamma = M_PI/4 + gamma_straight;
        }
        else if(gamma < - M_PI/4 + gamma_straight){
            gamma = - M_PI/4 + gamma_straight;
        } 

        //out << "phi : " << phi/M_PI << "*pi , gamma : " << gamma/M_PI <<"*pi"<< endl;
        delta_x = dist_*cos(phi)*sin(gamma);
        delta_y = dist_*sin(phi)*sin(gamma);
        delta_z = dist_*cos(gamma);

        //out << "delta_y : "<< delta_y << ", delta_x :"<< delta_x <<", delta_z :" << delta_z<< endl;
        
        pos_ = new_pos;
            
        pos_[0] +=  delta_x;
        pos_[1] +=  delta_y;
        pos_[2] +=  delta_z;

        // check distance with border

        if (!check_borders(pos_, rad*sqrt(1+ax.volume_inc_perc) + barrier_tickness)) {

            // check collision with other spheres
            Dynamic_Sphere s(pos_, rad,ax.volume_inc_perc,ax.swell, axon_id, 1, ax.active_state);

            if (!isSphereColliding(s, out)){  
                sphere_id += 4;
                prev_pos = new_pos;
                new_pos = pos_;
                centers.push_back(new_pos);
                achieved = true;
                sph_radii.push_back(rad);
                distances.push_back(dist_);
                return true;
                break;
            } 
            else {
                tries += 1;
            }
        }
        else {
            tries += 1;
         }
    }
    return false;
}

void AxonGammaDistribution::fiber_collapse(Eigen::Vector3d& new_pos, Eigen::Vector3d& prev_pos, std::vector<Eigen::Vector3d>& centers, int& fibre_collapsed_nbr, ostream& out){
    int n = centers.size();
    out << "Max tries reached ! centers.size() : " << centers.size() << endl;
    Eigen::Vector3d p = {centers[n-2-fibre_collapsed_nbr][0], centers[n-2-fibre_collapsed_nbr][1],centers[n-2-fibre_collapsed_nbr][2]};
    new_pos = p;
    out << "new_pos :" << p[2]<< endl;
    Eigen::Vector3d p2 ={centers[n-3-fibre_collapsed_nbr][0], centers[n-3-fibre_collapsed_nbr][1],centers[n-3-fibre_collapsed_nbr][2]};
    prev_pos = p2;
    out << "prev_pos :" <<p2[2] << endl;
    fibre_collapsed_nbr += 1;
    for (unsigned j=0; j< fibre_collapsed_nbr+1; ++j){
        centers.pop_back();
    }
}

void AxonGammaDistribution::fill_spheres_in_between(Axon ax, std::vector<Eigen::Vector3d>& centers, std::vector<Dynamic_Sphere>& spheres_to_add, std::vector<double>& sph_radii){

    Dynamic_Sphere last_sphere;
    for (unsigned i=0; i< centers.size(); ++i){
        if (i != 0){
            // squeleton has centers every r, so 3 spheres must be created in between

            Eigen::Vector3d center_between0 = {centers[i][0]/4+3*centers[i-1][0]/4, centers[i][1]/4+3*centers[i-1][1]/4, centers[i][2]/4+3*centers[i-1][2]/4};
            Dynamic_Sphere sphere0 (center_between0, sph_radii[i-1],ax.volume_inc_perc,  ax.swell, ax.id, 1, ax.active_state);  
            if (sphere0.center[2] < max_limits[2]){
                last_sphere = sphere0;
                spheres_to_add.push_back(sphere0);
            }

            Eigen::Vector3d center_between1 = {(centers[i][0]+centers[i-1][0])/2, (centers[i][1]+centers[i-1][1])/2, (centers[i][2]+centers[i-1][2])/2};
            Dynamic_Sphere sphere1(center_between1, sph_radii[i-1],ax.volume_inc_perc,  ax.swell, ax.id, 1, ax.active_state);  
            if (sphere1.center[2] < max_limits[2]){
                last_sphere = sphere1;
                spheres_to_add.push_back(sphere1);
            }

            Eigen::Vector3d center_between2 = {3*centers[i][0]/4+centers[i-1][0]/4, 3*centers[i][1]/4+centers[i-1][1]/4, 3*centers[i][2]/4+centers[i-1][2]/4};
            Dynamic_Sphere sphere2(center_between2, sph_radii[i-1],ax.volume_inc_perc,  ax.swell, ax.id, 1, ax.active_state);  
            if (sphere2.center[2] < max_limits[2]){
                last_sphere = sphere2;
                spheres_to_add.push_back(sphere2);
            }
        //out << sphere.center[0] << " " << sphere.center[1] << " " << sphere.center[2] << " " << sphere.radius << endl;
        }
        Dynamic_Sphere sphere(centers[i], sph_radii[i],ax.volume_inc_perc,  ax.swell, ax.id, 1, ax.active_state);  
        if (sphere.center[2] < max_limits[2]){
            last_sphere = sphere;
            spheres_to_add.push_back(sphere);
        }
        if(i ==centers.size()-1){

            //out << "last sphere : " << last_sphere.center[0] << " " << last_sphere.center[1] << " " << last_sphere.center[2] << " " << last_sphere.radius << endl;
            last_sphere.set_center( {last_sphere.center[0], last_sphere.center[1], max_limits[2]});
            spheres_to_add.push_back(last_sphere);
        }
    }
}

bool shrink_sphere_rad(double& rad, double& dist_){
    rad = rad *(1-0.01);
    dist_ = rad;
}

std::vector<Dynamic_Sphere> AxonGammaDistribution::GrowAxon(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out){

    std::vector<Eigen::Vector3d> centers;
    // add first position to the list of centers
    Eigen::Vector3d new_pos = ax.begin;
    centers.push_back(new_pos);
    int max_fibre_collapse = 50 ;
    std::vector<double> sph_radii;
    std::vector<double> distances;
    double dist_ = ax.radius;
    double rad = ax.radius;


    Eigen::Vector3d prev_pos = {new_pos[0], new_pos[1], new_pos[2]-ax.radius};

    Dynamic_Sphere s1(new_pos, ax.radius,ax.volume_inc_perc,ax.swell, axon_id, 1, ax.active_state);
    std::vector<Dynamic_Sphere> spheres_to_add;

    bool stop = false;
    int fibre_collapsed_nbr = 0;
    out << "try_axon :" << axon_id << endl;

    if(isSphereColliding(s1,out)) {
        out << "no room !" << endl;
        return spheres_to_add;
    }
    do{  

        bool next_center_found = find_next_center(ax, distances, dist_, sph_radii, rad, new_pos, prev_pos, centers,axon_id, out);

        // fibre collapse (see CONFIG)
        if(!next_center_found){
            
            if (centers.size() > 3+fibre_collapsed_nbr && fibre_collapsed_nbr < max_fibre_collapse){
     
                fiber_collapse(new_pos, prev_pos, centers, fibre_collapsed_nbr, out);
            }
            else{

                return spheres_to_add;
            }
        }

    }while (new_pos[2] < ax.end[2] + dist_);

    // create axon with 4 spheres in between each center
    fill_spheres_in_between(ax, centers, spheres_to_add, sph_radii);

    return spheres_to_add;
}   