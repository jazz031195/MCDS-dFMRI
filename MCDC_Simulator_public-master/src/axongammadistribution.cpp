#include "axongammadistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"
#include <chrono>


using namespace std;
using namespace Eigen;
using namespace std::chrono;

AxonGammaDistribution::AxonGammaDistribution (double dyn_perc_,double volume_inc_perc_, unsigned num_ax, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l, float min_radius, bool active_state_, double c2_)
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
    duration = 0.0;
    c2 = c2_;
    tortuosities.clear();


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

void AxonGammaDistribution::find_target_point (double c2, double radius, Eigen::Vector3d& initial_point , Eigen::Vector3d& target_point){

    bool achieved = false;

    while(!achieved){
        std::random_device rd;
        std::mt19937 gen(rd()); 
        std::normal_distribution<float> cos_2_dist (c2, 0.04); 

        double cos_2 = cos_2_dist(gen);
        while (cos_2 < 0 || cos_2 > 1){
            cos_2 = cos_2_dist(gen);
        }
        double sin_2 = 1-cos_2;
        double slope_2 = sin_2/cos_2;
        double diff_square_x_y = slope_2*max_limits[2];

        std::uniform_real_distribution<double> udist(0, diff_square_x_y);
        double diff_y_2 = udist(gen);
        double diff_y = sqrt(diff_y_2);
        double diff_x_2 = diff_square_x_y - diff_y_2;
        double diff_x = sqrt(diff_x_2);
        // 0: x,y > 0, 1: x > 0 and y < 0, 2: x<0 and y>0, 3 : x,y <0
        std::uniform_real_distribution<float> pos_or_neg (0, 4);
        int pos_or_neg_ = int(pos_or_neg(gen));
        if (pos_or_neg_ == 1 || pos_or_neg_== 3){
            diff_y = -diff_y;
        }
        if(pos_or_neg_ == 2 || pos_or_neg_== 3){
            diff_x = -diff_x;
        }

        std::uniform_real_distribution<double> udist_(0, 1);
        double t = udist_(gen);
        double distance_to_border = radius*sqrt(1+volume_inc_perc) + barrier_tickness;
        double x_i = (t * (max_limits[0]-distance_to_border)) + (1 - t) * (min_limits[0]+ distance_to_border);
        t = udist_(gen);
        double y_i = (t * (max_limits[1]-distance_to_border) + (1 - t) * (min_limits[1]+ distance_to_border));

        double x_f = x_i + diff_x;
        double y_f = y_i + diff_y;

        if (check_borders({x_f, y_f, max_limits[2]}, radius*sqrt(1+volume_inc_perc))){
            achieved = true;
            initial_point = {x_i, y_i, min_limits[2]};
            target_point = {x_f, y_f, max_limits[2]};
        }
    }
}


void AxonGammaDistribution::createGammaSubstrate(ostream& out)
{
    // generate the gamma distribution
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    uint repetition = 1;
    uint max_adjustments = 0;
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

    auto start = high_resolution_clock::now();

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


    while (!achieved)
    {

        double target_icvf = this->icvf - adjustments * adj_increase;
        // computes max_limits
        computeMinimalSize(radiis, target_icvf, max_limits);


        for (uint t = 0; t < repetition; t++)
        {
            //vector<Axon> axons_to_add;

            axons.clear();

            unsigned stuck = 0;

            Vector3d Q;
            Vector3d D;

            for (unsigned i = 0; i < num_obstacles; i++)
            {
                if(stuck < 3000){

                    if (c2 == 1){
                    double t = udist(gen);
                    double distance_to_border = radiis[i]*sqrt(1+volume_inc_perc) + barrier_tickness;
                    double x = (t * (max_limits[0]-distance_to_border)) + (1 - t) * (min_limits[0]+ distance_to_border);
                    t = udist(gen);
                    double y = (t * (max_limits[1]-distance_to_border) + (1 - t) * (min_limits[1]+ distance_to_border));

                    //out << " obstacle :"<< i << endl;

                    Q = {x, y, min_limits[2]};
                    D = {x, y, max_limits[2]};
                    }
                    else{

                        find_target_point (c2, radiis[i], Q, D);

                    }

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
                break;
            }
        }
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
    icvf_current = best_icvf;

    auto stop = high_resolution_clock::now();
    auto duration_ = duration_cast<seconds>(stop - start);
    duration = duration_.count() ;


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
        out << "Axon: " << i << " tortuosity:"<< tortuosities[i] << endl;
    }
    out << "Time_to_grow:" << duration  << "_seconds"<< endl;
}

double AxonGammaDistribution::computeICVF(std::vector<Axon> &axons, Vector3d &min_limits, Vector3d &max_limits, int &num_no_repeat)
{

    if (axons.size() == 0)
        return 0;

    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1])*(max_limits[2] - min_limits[2]);

    double AreaC = 0;

    double tortuosity;

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

        double ax_length = 0;

        double mean_rad= 0;
        double rads  =0 ;

        if (axons[i].spheres.size() > 0){

            for (uint j = 0; j < axons[i].spheres.size(); j++){
                if (j > 0){
                    double l = (axons[i].spheres[j-1].center-axons[i].spheres[j].center).norm();
                    ax_length += l;
                }
                mean_rad += axons[i].spheres[j].radius;

            }

            mean_rad = mean_rad/axons[i].spheres.size();

            tortuosity = ax_length/((axons[i].begin-axons[i].end).norm());

            tortuosities.push_back(tortuosity);

            AreaC += M_PI * mean_rad * mean_rad * ax_length;
        }
        else{
            double rad = axons[i].radius;
            AreaC += M_PI *rad * rad * ax_length;
        }
        
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

    //Eigen::Vector3d straight_vector = {new_pos[0] - prev_pos[0], new_pos[1] - prev_pos[1], new_pos[2] - prev_pos[2]};
    //Eigen::Vector2d straight_vector_xy = {straight_vector[0] ,straight_vector[1] };
    //straight_vector_xy = straight_vector_xy.normalized();
    //Eigen::Vector2d straight_vector_xz = {straight_vector[0] ,straight_vector[2] }; 
    //straight_vector_xz = straight_vector_xz.normalized();

    //out << "straight_vector_xy : ( " << straight_vector_xy[0] << ", " << straight_vector_xy[1] << " )" << endl;
    //out << "straight_vector_xz : ( " << straight_vector_xz[0] << ", " << straight_vector_xz[1] << " )" << endl;


    double phi_to_target = atan2(vector_to_target_xy[1], vector_to_target_xy[0]);
    double gamma_to_target = atan2(vector_to_target_xz[0],vector_to_target_xz[1]);

    // caculate angle wrt x axis
    //double phi_straight = atan2(straight_vector_xy[1],straight_vector_xy[0]);
    // calculate angle wrt z axis
    //double gamma_straight = atan2(straight_vector_xz[0],straight_vector_xz[1]);

    return std::make_tuple(phi_to_target, gamma_to_target);
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

bool AxonGammaDistribution::isSphereColliding(Dynamic_Sphere sph){
    
    Vector3d position = sph.center;
    double distance_to_be_inside = 10* barrier_tickness + sph.max_radius;
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

bool AxonGammaDistribution::find_next_center(Axon ax, Dynamic_Sphere& s, double dist_, double& rad, Eigen::Vector3d& new_pos, Eigen::Vector3d& prev_pos, int axon_id, ostream& out){

    double phi;
    double gamma;
    double phi_to_target, gamma_to_target;
    bool achieved = false;
    int tries = 0;
    int max_tries = 50000;
    int sphere_id = 0;
    Eigen::Vector3d pos_;

    std::random_device rd;
    std::mt19937 gen(rd()); 

    double delta_x;
    double delta_y;
    double delta_z;

    double angle_between = M_PI/3;

    while(!achieved && tries < max_tries){
        tie(phi_to_target, gamma_to_target) = phi_gamma_to_target (prev_pos, new_pos, ax.end,  out);

        while (angle_between >= M_PI/3){
            std::normal_distribution<float> phi_dist (phi_to_target/M_PI, 0.1); 
            phi = phi_dist(gen)*M_PI;
            //phi = phi_to_target;
            std::normal_distribution<float> gamma_dist (gamma_to_target/M_PI, 0.1); 
            gamma = gamma_dist(gen)*M_PI;
            //gamma = gamma_to_target;

            //out << "phi : " << phi/M_PI << "*pi , gamma : " << gamma/M_PI <<"*pi"<< endl;
            delta_x = dist_*cos(phi)*sin(gamma);
            delta_y = dist_*sin(phi)*sin(gamma);
            delta_z = dist_*cos(gamma);

            //out << "delta_y : "<< delta_y << ", delta_x :"<< delta_x <<", delta_z :" << delta_z<< endl;
            
            pos_ = new_pos;
            Eigen::Vector3d prev_to_pos = (new_pos-prev_pos).normalized();
                
            pos_[0] +=  delta_x;
            pos_[1] +=  delta_y;
            pos_[2] +=  delta_z;

            Eigen::Vector3d pos_to_new_pos = (pos_-new_pos).normalized();

            angle_between = acos(pos_to_new_pos.dot(prev_to_pos));

            if(angle_between >= M_PI/3){
                out << "angle: " << angle_between/M_PI <<  endl;
            }

            //out << "angle between :" << angle_between/M_PI << " PI " << endl;
        }

        // check distance with border
        if (!check_borders(pos_, rad*sqrt(1+ax.volume_inc_perc) + barrier_tickness)) {

            // check collision with other spheres
            Dynamic_Sphere sphere_ (pos_, rad,ax.volume_inc_perc,ax.swell, axon_id, 1, ax.active_state);

            if (!isSphereColliding(sphere_)){ 
                s =sphere_; 
                achieved = true;
                //out << "sphere at center : "<< sphere_.center<< ",with radius "<< sphere_.min_radius <<" doesn't collide" << endl;
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
    int nbr_discarded_centers = (fibre_collapsed_nbr+1)*4;
    Eigen::Vector3d p = {centers[n-nbr_discarded_centers-1][0], centers[n-nbr_discarded_centers-1][1],centers[n-nbr_discarded_centers-1][2]};
    new_pos = p;
    Eigen::Vector3d p2 ={centers[n-nbr_discarded_centers-5][0], centers[n-nbr_discarded_centers-5][1],centers[n-nbr_discarded_centers-5][2]};
    prev_pos = p2;
    for (unsigned j=0; j< nbr_discarded_centers; ++j){
        centers.pop_back();
    }
    fibre_collapsed_nbr += 1;
}

bool AxonGammaDistribution::fill_spheres_in_between(Axon ax, Dynamic_Sphere& added_sphere, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii, ostream& out){
    std::vector<Dynamic_Sphere> spheres_to_add_;
    int centers_size = centers.size();
    Eigen::Vector3d last_center = {centers[centers_size-1][0], centers[centers_size-1][1],centers[centers_size-1][2]};
    double last_radius = sph_radii[centers_size-1];
    

    for (unsigned i = 0; i < 4; i++){
        Eigen::Vector3d center_between = {(3-i)*last_center[0]/4+(i+1)*added_sphere.center[0]/4, (3-i)*last_center[1]/4+(i+1)*added_sphere.center[1]/4, (3-i)*last_center[2]/4+(i+1)*added_sphere.center[2]/4};
        double radius_in_between = (3-i)*last_radius/4+(i+1)*added_sphere.min_radius/4;
        Dynamic_Sphere sphere (center_between, radius_in_between,ax.volume_inc_perc, ax.swell, ax.id, 1, ax.active_state);  
        if (i== 3){
            spheres_to_add_.push_back(sphere);
        }
        else if (!isSphereColliding(sphere)){
            //out << "overlapping sphere :" << i << " is ok " << endl;
            spheres_to_add_.push_back(sphere);
        }

    }
    if(spheres_to_add_.size() == 4){
        for (unsigned i = 0; i < spheres_to_add_.size(); i++){
            centers.push_back(spheres_to_add_[i].center);
            sph_radii.push_back(spheres_to_add_[i].min_radius);
        }

        return true;
    }
    else{
        return false;
    }
}

void AxonGammaDistribution::shrink_sphere_rad(double& rad, double axon_rad, double& shrink_perc, ostream& out){
    shrink_perc += 0.02;
    rad = axon_rad *sqrt(1-shrink_perc);
}



std::vector<Dynamic_Sphere> AxonGammaDistribution::GrowAxon(Axon ax, double distance_to_be_inside, int axon_id,  ostream& out){

    std::vector<Eigen::Vector3d> centers;
    // add first position to the list of centers
    Eigen::Vector3d new_pos = ax.begin;
    int max_fibre_collapse = 5;
    double max_shrinking = 0.7;
    std::vector<double> sph_radii;
    double rad = ax.radius;
    // initial distance between spheres is the radius
    double dist_ = rad;
    double prev_radius;


    Eigen::Vector3d prev_pos = {new_pos[0], new_pos[1], new_pos[2]-dist_};

    Dynamic_Sphere s1(new_pos, ax.min_radius,ax.volume_inc_perc,ax.swell, axon_id, 1, ax.active_state);
    std::vector<Dynamic_Sphere> spheres_to_add;

    bool stop_fibre_collapsing = false;
    int fibre_collapsed_nbr = 0;
    double shrink_perc = 0.0;

    out << "try_axon :" << axon_id << endl;

    if(isSphereColliding(s1)) {
        return spheres_to_add;
    }
    else{
        centers.push_back({s1.center[0],s1.center[1],s1.center[2]});
        sph_radii.push_back(s1.min_radius);
        prev_radius = s1.min_radius;
    }
    do{
        Dynamic_Sphere sphere_to_add;

        bool sphere_filling = false;
        // find a sphere to add that doesn't collide with environment at dist_ specified, add to sphere_to_add
        bool next_center_found = find_next_center(ax, sphere_to_add, dist_, rad, new_pos, prev_pos, axon_id, out);
        if (next_center_found){
            // checks if all spheres in between that overlap do not collide with environment, add to centers and sph_radii
            sphere_filling = fill_spheres_in_between(ax, sphere_to_add, centers, sph_radii, out);
        }

        if(!sphere_filling){
            if (!stop_fibre_collapsing && (centers.size() > (5+(fibre_collapsed_nbr+1)*4)) && (fibre_collapsed_nbr < max_fibre_collapse)){
                // fibre collapse (see CONFIG)
                //out << "fiber collapse nbr :" << fibre_collapsed_nbr << endl;
                fiber_collapse(new_pos, prev_pos, centers, fibre_collapsed_nbr, out);
            }
            else{
            // if another sphere cannot be added, we want to check if by shrinking it at the maximum value it can pass

                stop_fibre_collapsing = true;
            // if a very small sphere can be added, we want to find the largest one that can be added
                if(rad > min_radius && shrink_perc + 0.01 < max_shrinking ){
                    // changes rad and dist_ values
                    // rad is shrinked and dist_ is recalculated
                    shrink_sphere_rad(rad, ax.radius, shrink_perc, out);
                    dist_ = max(rad, prev_radius);
                    //out << "shrink radius by " << shrink_perc*100 << " percent" << endl;
                }
                else{

                    return spheres_to_add;
                }
            }

        }
        else{

            // reinitialise radius and distance between spheres
            
            prev_radius = rad;
            dist_ = ax.radius;  
            rad = ax.radius;
            shrink_perc = 0.0;
            new_pos = {centers[centers.size() -1][0],centers[centers.size() -1][1],centers[centers.size() -1][2]};
            prev_pos = {centers[centers.size()-3][0],centers[centers.size() -3][1],centers[centers.size() -3][2]};
        }
        //out << "centers.size() :" << centers.size() <<endl;

    //out << "new_pos[2] : " << new_pos[2] <<", ax.end[2] :" << ax.end[2] << endl;
    }while (new_pos[2] < ax.end[2]);

    // add spheres to spheres_to_add
    out << "spheres :" << endl;
    bool last_sphere_add =false;
    for (unsigned i=0; i< sph_radii.size(); ++i){
        if (!last_sphere_add && centers[i][2]<= max_limits[2]){
            Dynamic_Sphere s(centers[i], sph_radii[i],ax.volume_inc_perc,ax.swell, axon_id, 1, ax.active_state);
            spheres_to_add.push_back(s);
            out << "sphere at : " << centers[i] << endl;
        }
        else if (!last_sphere_add){
            double dist_i_1= abs((centers[i-1] - max_limits).norm());
            double dist_i = abs((centers[i] - max_limits).norm());
            Eigen::Vector3d vector = (centers[i] - centers[i-1]).normalized();
            // centers[i-1][2] + lambda*vector[2] = max_limits[2]
            double lambda_ = (max_limits[2]-centers[i-1][2])/vector[2];
            Eigen::Vector3d last_center = centers[i-1] +lambda_*vector;
            double weight_i_1 = abs((last_center-centers[i]).norm())/abs((centers[i-1]-centers[i]).norm());
            double weight_i = abs((last_center-centers[i-1]).norm())/abs((centers[i-1]-centers[i]).norm());
            double last_radius = weight_i*sph_radii[i]+weight_i_1*sph_radii[i-1];
            Dynamic_Sphere last_sphere(last_center, last_radius,ax.volume_inc_perc,ax.swell, axon_id, 1, ax.active_state);
            last_sphere_add = true;
            if (!isSphereColliding(last_sphere)){
                spheres_to_add.push_back(last_sphere);
                out << "last sphere add : " << last_center << endl;
            }
        }
        if(last_sphere_add){
            out << "distance to target :" << (centers[centers.size()-1]-ax.end).norm() << endl;
            return spheres_to_add;
        }
        //out << " radius " << i << " :" << sph_radii[i] << endl;
    }


    return spheres_to_add;

}   