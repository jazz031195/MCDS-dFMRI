#include "axongammadistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"
#include <chrono>



using namespace std;
using namespace Eigen;
using namespace std::chrono;

AxonGammaDistribution::AxonGammaDistribution (double dyn_perc_,double volume_inc_perc_, unsigned num_ax, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l, double min_radius, bool active_state_, double c2_, bool tortuous_, double step_length_)
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
    tortuous = tortuous_;
    small_voxel_size = max_limits;
    step_length = step_length_;

}
void AxonGammaDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d &l)
{
    /*
    Computes the minimal Voxel size for the chosen icvf, by assuming straight axons

    radiis : list of radii for each axon to grow
    icvf_ : Intracompartment volume fraction
    l : 3D vector with each of the 3 values being the length of the vector
    */

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
    /* 
    Displays the Gamma Distribution of the axons
    */
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
        SimErrno::info(message,std::cout);
    }
    message = "9-10:" + std::string(p[9]*nstars/nrolls,'*') ;
    SimErrno::info(message,std::cout);
    message = ">10: " +  std::string(p[10]*nstars/nrolls,'*') + "\n" ;
    SimErrno::info(message,std::cout);
}

void AxonGammaDistribution::find_target_point (double c2, double radius, Eigen::Vector3d& initial_point , Eigen::Vector3d& target_point){
    /*
    Finds the target of the axon (position of last sphere) w.r.t the c2 value given.
    
    c2 : Mean(cos²(angle))
    radius : Radius of axon
    initial_point : Position from which the axon will start to grow
    target_point : Position of the target of the growth
    */

    // find initial point in random position in voxel
    std::random_device rd;
    std::mt19937 gen(rd()); 
    bool achieved = false;
    std::uniform_real_distribution<double> udist_(0, 1);
    double t = udist_(gen);
    //double distance_to_border = radius*sqrt(1+volume_inc_perc) + barrier_tickness;
    double x_i = (t * (small_voxel_size[0])) + (1 - t) * (min_limits[0]);
    t = udist_(gen);
    double y_i = (t * (small_voxel_size[1]) + (1 - t) * (min_limits[1]));


    while(!achieved){

        // c2 is mean(cos²(angle)), so we take cos²(angle) from a distribution that has c2 as mean
        std::normal_distribution<float> cos_2_dist (0, 1-c2); 
        double cos_2 = cos_2_dist(gen);
        while (cos_2 < 0 || cos_2 > 1){
            cos_2 = cos_2_dist(gen);
        }

        // target_point - initial_point = diff
        // diff_2 = diff²
        double diff_x_2, diff_y_2, diff_z_2;
        // vi.dot(diff) = cos(angle) (between vo (straight vector) and diff (vector in angle))
        // vi = [0,0,1] -> diff[z] = cos(angle)
        // cos²(angle) = cos_2 = diff_z²
        diff_z_2 = cos_2;
        std::uniform_real_distribution<double> udist_y(diff_z_2-1, 1-diff_z_2);
        diff_y_2 = udist_y(gen);
        std::uniform_real_distribution<double> udist_x(diff_z_2+diff_y_2 -1, 1-diff_z_2-diff_y_2);
        diff_x_2 = udist_x(gen);

        double length = max_limits[2]/sqrt(diff_z_2);
        double x_f = x_i + length*sqrt(diff_x_2);
        double y_f = y_i + length*sqrt(diff_y_2);

        // check if target point is in the voxel limits 
        Eigen::Vector2d twin_delta_pos = {0,0};
        if (!check_borders({x_f, y_f, max_limits[2]}, radius*sqrt(1+volume_inc_perc)+step_length, twin_delta_pos)){
            achieved = true;
            initial_point = {x_i, y_i, min_limits[2]};
            target_point = {x_f, y_f, max_limits[2]};
        }
    }
}

void display_progress(double nbr_axons, double number_obstacles){
    int cTotalLength = 50;
    double lProgress = nbr_axons/number_obstacles;
    std::cout << 
        "\r[" <<                                            //'\r' aka carriage return should move printer's cursor back at the beginning of the current line
            string(int(cTotalLength * lProgress), '*') <<        //printing filled part
            string(int(cTotalLength * (1 - lProgress)), '-') <<  //printing empty part
        "] "  << nbr_axons << "/" << number_obstacles << endl;          //printing percentage
} 

void AxonGammaDistribution::flip(int flip_nbr, int j, int k, Eigen::Vector3d small_voxel_size, Eigen::Vector3d& initial_pos){

    // projection onto next small voxel
    double x_increment = j*small_voxel_size[0];
    double y_increment = k*small_voxel_size[1];
    std::cout << "   small_voxel_size: " << small_voxel_size << endl;
    std::cout << "initial position : " << initial_pos << endl;
    initial_pos = {initial_pos[0] + x_increment, initial_pos[1] + y_increment, initial_pos[2]};
    std::cout << " position after increment : " << initial_pos << endl;
    double x_flip;
    double y_flip;
    // center of small voxel
    Eigen::Vector2d center = {(j+0.5)*small_voxel_size[0], (k+0.5)*small_voxel_size[1]};
    std::cout << "center : " << center << endl;

    std::cout << "flip nbr : " << flip_nbr << endl;
    if (flip_nbr == 0){
        x_flip = 0;
        y_flip = 0;
    }  
    // flips of the small voxel wrt y
    else if (flip_nbr == 1){
        double distance = 2*(center[1]-initial_pos[1]);
        y_flip = distance;
    }
    // flips of the small voxel wrt x
    else if (flip_nbr == 2){
        double distance = 2*(center[0]-initial_pos[0]);
        x_flip = distance;
    }
    else if (flip_nbr == 3){
        double distance = 2*(center[0]-initial_pos[0]);
        x_flip = distance;
        distance = 2*(center[1]-initial_pos[1]);
        y_flip = distance;
    } 

    initial_pos = {initial_pos[0] + x_flip, initial_pos[1] + y_flip, initial_pos[2]};
    std::cout << " position after flip : " << initial_pos << endl;
}

void AxonGammaDistribution::add_periodic_voxel(int nbr_small_voxels, Eigen::Vector3d small_voxel_size){
    
    std::vector<Axon> axons_to_add;
    int flip_nbr;
    for (unsigned j = 0; j < nbr_small_voxels -1 ; ++j){
        
        for (unsigned k = 0; k < nbr_small_voxels -1 ; ++k){
            
            int flip_nbr = 0; 
            
            for (unsigned i = 0; i < axons.size(); ++i){
                if (j == 0 and k == 0){
                    k += 1;
                }
                Eigen::Vector3d new_begin = axons[i].begin;
                flip(flip_nbr, j, k, small_voxel_size, new_begin);
                Eigen::Vector3d new_end = axons[i].end;
                flip(flip_nbr, j, k, small_voxel_size, new_end);
                Axon *ax = new Axon (axons[i]);
                ax->begin = new_begin;
                ax->end = new_end;
                std::vector<Dynamic_Sphere> new_spheres;
                bool stop = false;
                for (unsigned s = 0; s < axons[i].spheres.size(); ++s){
                    if (!stop){ 
                        Dynamic_Sphere new_sphere = axons[i].spheres[s];
                        // if axon already exists
                        if (isSphereColliding(new_sphere)){
                            stop = true;
                        } 
                        Eigen::Vector3d new_center = new_sphere.center;
                        flip(flip_nbr, j, k, small_voxel_size, new_center);
                        new_sphere.center = new_center;
                        new_spheres.push_back(new_sphere);
                    } 
                }
                if (new_spheres.size() != 0){  
                    ax->set_spheres(new_spheres);
                    axons_to_add.push_back(*ax);
                }
                else{
                    delete ax;
                }  
            }
        }
    }
    // extend axons with axons to add
    axons.reserve(axons.size() + distance(axons_to_add.begin(),axons_to_add.end()));
    axons.insert(axons.end(),axons_to_add.begin(),axons_to_add.end());
}

void AxonGammaDistribution::createGammaSubstrate(ostream& out)
{
    /* 
        Generates the gamma distribution of axons. 
    */
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    uint repetition = 1;
    min_limits = {0., 0., 0.};

    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0, 1);
    std::vector<double> radiis(num_obstacles, 0);

    int number_swelling_axons = int(num_obstacles * dyn_perc);
    std::vector<bool> bool_swell_ax_id(num_obstacles, false);

    bool achieved = false;

    int tried = 0;

    string message;

    auto start = high_resolution_clock::now();
    std::cout << " Growing axons " << endl;

    // create list bools that show whether an axon has the potential to swell
    for (unsigned i = 0; i < number_swelling_axons; ++i)
    {
        int random_id = rand() % num_obstacles;
        if (i > 0)
        {
            while (bool_swell_ax_id[random_id]) {
                random_id = rand() % num_obstacles;
            }
        }
        bool_swell_ax_id[random_id] = true;
    }

    for (unsigned i = 0; i < num_obstacles; ++i)
    {
        
        if (tried > 1000)
        {
            message = " Radii distribution cannot be sampled [Min. radius Error]\n";
            SimErrno::error(message, std::cout);
            assert(0);
        }
        double jkr = distribution(generator);

        // generates the radii in a list
        if (jkr < 0.05)
        {
            i--;
            tried++;
            continue;
        }
        tried = 0;

        radiis[i] = jkr * 1e-3; // WE CONVERT FROM UM TO MM HERE
    }

    // sorts the radii 
    std::sort(radiis.begin(), radiis.end(), [](const double a, double b) -> bool
              { return a > b; });

    double max_radius_ = radiis[0] * sqrt(1+volume_inc_perc);

    int number_axons_placed = 0;

    while (!achieved)
    {

        computeMinimalSize(radiis, icvf, small_voxel_size);
        // in case the voxel size input in parameters is too small 
        if (small_voxel_size[2] > max_limits[2]){
            max_limits = small_voxel_size;
        }
        // count how many times this small voxel can fit in the voxel size initialised
        int sqrt_nbr_small_voxels = int((max_limits[0]/small_voxel_size[0]));
        // readjust max_limits so that it is a multiple of small_voxel
        max_limits = {small_voxel_size[0]*sqrt_nbr_small_voxels, small_voxel_size[1]*sqrt_nbr_small_voxels, small_voxel_size[2]*sqrt_nbr_small_voxels};
        // make this small voxel a rectangle
        small_voxel_size[2] = max_limits[2];
        axons.clear();

        unsigned stuck = 0;

        Vector3d Q;
        Vector3d D;

        int number_twins = 0;
        

        for (unsigned i = 0; i < num_obstacles; i++)
        {
            bool axon_placed = false;
            if(stuck < 1000000 && !axon_placed){

                if (c2 == 1){
                    double t = udist(gen);
                    //double distance_to_border = radiis[i]*sqrt(1+volume_inc_perc) + step_length;
                    double x = (t * (small_voxel_size[0])) + (1 - t) * (min_limits[0]);
                    t = udist(gen);
                    double y = (t * (small_voxel_size[1]) + (1 - t) * (min_limits[1]));

                    Q = {x, y, min_limits[2]};
                    D = {x, y, small_voxel_size[2]};
                }
                else{

                    find_target_point (c2, radiis[i], Q, D);

                }

                Axon *ax = new Axon (number_axons_placed, radiis[i], Q, D, volume_inc_perc, active_state, bool_swell_ax_id[i], 1);
                out << "axon : " << ax->id << " , i :" << i << endl;
                bool is_outside_subvoxel = false;
                std::vector<Eigen::Vector2d> all_twin_delta_pos;
                bool can_shrink = false;
                if(stuck > 100) {
                    can_shrink = true;
                }
                std::vector<Dynamic_Sphere> spheres_to_add = GrowAxon(ax, max_radius_, is_outside_subvoxel, all_twin_delta_pos, can_shrink, out);
                if(spheres_to_add.size() != 0)
                {
                    ax->set_spheres(spheres_to_add);
                    axons.push_back(*ax);
                    number_axons_placed += 1;
                    if(is_outside_subvoxel) {
                        vector<Axon*> twin_axons;
                        out << "all_twin_delta_pos.size() :" << all_twin_delta_pos.size() << endl;
                        for (unsigned k = 0; k < all_twin_delta_pos.size(); k++){
                            Eigen::Vector2d twin_delta_pos = all_twin_delta_pos[k];
                            out << "twin_delta_pos :" << twin_delta_pos << endl;
                            twin_axons.clear();
                            createTwinAxons(ax, twin_delta_pos, twin_axons, number_axons_placed);
                            out << "twin_axons :" << twin_axons.size() << endl;
                            for (unsigned j = 0; j < twin_axons.size(); j++){
                                out << "adding twin axon :" << twin_axons[j]->id << endl;
                                axons.push_back(*twin_axons[j]);
                                number_axons_placed += 1;
                                number_twins += 1;
                            }
                        }
                        
                    }
                    stuck = 0;
                    axon_placed = true;
                    
                }
                else{
                    i -= 1; 
                    stuck += 1;
                    delete ax;
                    continue;
                }
                display_progress(axons.size(), num_obstacles+number_twins);
            }


        } // end for axons

        // check ICVF
        std::cout << "sqrt_nbr_small_voxels :" << sqrt_nbr_small_voxels << endl;
        add_periodic_voxel(sqrt_nbr_small_voxels, small_voxel_size);
        icvf_current = computeICVF();
        achieved = true;
    }
    num_obstacles = number_axons_placed;



    for (unsigned i = 0; i < axons.size(); i++){
        out << "Axon :" << i << endl;
        for (unsigned s = 0; s < axons[i].spheres.size(); s++){
            out << axons[i].spheres[s].center[0] << " " << axons[i].spheres[s].center[1] << " " << axons[i].spheres[s].center[2] << " " << axons[i].spheres[s].radius << endl;
        }
    }

    // time
    auto stop = high_resolution_clock::now();
    auto duration_ = duration_cast<seconds>(stop - start);
    duration = duration_.count() ;

    // messages
    out <<"icvf:"<< icvf_current << "voxel size: "<< max_limits[0] << endl;
    message = "ICVF achieved: " + to_string(icvf_current * 100) + "  (" + to_string(int((icvf_current / icvf * 100))) + "% of the desired icvf)\n";
    SimErrno::info(message, std::cout);

}

void AxonGammaDistribution::printSubstrate(ostream &out)
{
    double scale = 1;
    out << scale << endl;
    out << volume_inc_perc << endl;
    out << dyn_perc << endl;
    out << icvf_current << endl;
    out << min_limits[2] << endl;
    out << max_limits[2] << endl;
    
    for (unsigned i = 0; i < axons.size(); i++)
    {
        for (unsigned s = 0; s < axons[i].spheres.size(); s++){

            out<< axons[i].spheres[s].center[0]/scale << " " << axons[i].spheres[s].center[1]/scale << " " << axons[i].spheres[s].center[2]/scale<< " " << axons[i].spheres[s].min_radius/scale << " "<<axons[i].spheres[s].swell << endl;
        }
        out << "Axon: " << axons[i].id << " tortuosity:"<< tortuosities[i] << endl;
    }
    out << "Time_to_grow:" << duration  << "_seconds"<< endl;
}

double AxonGammaDistribution::computeICVF()
{

    if (axons.size() == 0)
        return 0;

    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1])*(max_limits[2] - min_limits[2]);
    std::cout <<"volume tot : " << AreaV << endl;
    double AreaC = 0;

    double tortuosity;

    // using a lambda function:
    std::sort(axons.begin(), axons.end(), [](const Axon a, Axon b) -> bool
              { return a.radius > b.radius; });


    for (uint i = 0; i < axons.size(); i++)
    {

        double ax_length = 0;
        double mean_rad= 0;
        int num_spheres = 0;

        // if twin 
        if (i > 0 && axons[i].radius ==  axons[i-1].radius){
            tortuosities.push_back(tortuosities[-1]);
            continue;
        } 
         

        if (axons[i].spheres.size() > 1){

            for (uint j = 1; j < axons[i].spheres.size(); j++){           
                double l = (axons[i].spheres[j-1].center-axons[i].spheres[j].center).norm();
                ax_length += l;
                mean_rad += axons[i].spheres[j].min_radius;
                num_spheres += 1;

            }
            if (num_spheres != 0) { 
                mean_rad = mean_rad/num_spheres;
            }
            else{
                mean_rad = 0;
            }

            tortuosity = ax_length/((axons[i].begin-axons[i].end).norm());

            tortuosities.push_back(tortuosity);

            AreaC += M_PI * mean_rad * mean_rad * ax_length;
        }
        
    }
    return AreaC / AreaV;
}



std::tuple<double, double>  phi_theta_to_target (Eigen::Vector3d new_pos, Eigen::Vector3d end,  ostream& out) {
    /*
    Find phi and theta angles in spherical coordinates
    see : http://douglashalse.com/index.php/2019/09/25/spherical-coordinates/

    */
    Eigen::Vector3d vector_to_target = {end[0] - new_pos[0], end[1] - new_pos[1], end[2]+0.01 - new_pos[2]};
    vector_to_target = vector_to_target.normalized();
    double phi_to_target;
    double theta_to_target;

    if (vector_to_target[0] != 0){
        theta_to_target = atan(vector_to_target[1]/vector_to_target[0]);
    }
    else if (vector_to_target[0] == 0 && vector_to_target[1] == 0){
        theta_to_target = 0;
    }
    else if (vector_to_target[0] == 0){
        theta_to_target = M_PI/2;
    }

    if (vector_to_target[2] != 0){
        phi_to_target = atan(sqrt(vector_to_target[0]*vector_to_target[0]+vector_to_target[1]*vector_to_target[1])/vector_to_target[2]);
    }
    else{
        phi_to_target = M_PI/2;
    }
    return std::make_tuple(phi_to_target, theta_to_target);
}

bool AxonGammaDistribution::check_borders(Eigen::Vector3d pos, double distance_to_border, Eigen::Vector2d& twin_delta_pos) {
    
    
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
bool AxonGammaDistribution::checkAxoninBorders(Axon* ax, Eigen::Vector2d& twin_delta_pos){
    
    for (unsigned i = 0; i < ax->spheres.size() ; i++){
        if (check_borders(ax->spheres[i].center, ax->spheres[i].max_radius, twin_delta_pos) ){
            return false;
            break;
        } 
    } 
    return true;
} 

bool AxonGammaDistribution::checkAxonsinBorders(std::vector<Axon*> new_axons, Eigen::Vector2d& twin_delta_pos){
    
    for (unsigned i = 0; i < new_axons.size() ; i++){
        if (!checkAxoninBorders(new_axons[i], twin_delta_pos)){
            return false;
            break;
        } 
    } 
    return true;
} 

void AxonGammaDistribution::createTwinSpheres(std::vector<Dynamic_Sphere>& twin_spheres, Dynamic_Sphere s, Eigen::Vector2d twin_delta_pos){

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
  

void AxonGammaDistribution::createTwinAxons(Axon* ax, Eigen::Vector2d twin_delta_pos, std::vector<Axon*>& twin_axons, int id){
    
    Axon* twin_ax = new Axon (id, ax->min_radius, ax->begin, ax->end, ax->volume_inc_perc, ax->active_state, ax->swell, 1);
    std::vector<Dynamic_Sphere> twin_spheres;

    for (unsigned i = 0; i < ax->spheres.size() ; i++){
        Dynamic_Sphere twin_sphere = ax->spheres[i];
        twin_sphere.center = {twin_sphere.center[0]+twin_delta_pos[0], twin_sphere.center[1]+twin_delta_pos[1], twin_sphere.center[2]};
        twin_spheres.push_back(twin_sphere);
    }
    twin_ax->set_spheres(twin_spheres);
    twin_axons.push_back(twin_ax);

    // if axon overlaps with 2 different planes of boundary
    if (twin_delta_pos[0] != 0 && twin_delta_pos[1] != 0){
        Axon* twin_ax1 = new Axon (id+1, ax->min_radius, ax->begin, ax->end, ax->volume_inc_perc, ax->active_state, ax->swell, 1);
        twin_spheres.clear();

        for (unsigned i = 0; i < ax->spheres.size() ; i++){
            Dynamic_Sphere twin_sphere = ax->spheres[i];
            twin_sphere.center = {twin_sphere.center[0], twin_sphere.center[1]+twin_delta_pos[1], twin_sphere.center[2]};
            twin_spheres.push_back(twin_sphere);
        }
        twin_ax1->set_spheres(twin_spheres);
        twin_axons.push_back(twin_ax1);

        Axon* twin_ax2 = new Axon (id+2, ax->min_radius, ax->begin, ax->end, ax->volume_inc_perc, ax->active_state, ax->swell, 1);
        twin_spheres.clear();
        for (unsigned i = 0; i < ax->spheres.size() ; i++){
            Dynamic_Sphere twin_sphere = ax->spheres[i];
            twin_sphere.center = {twin_sphere.center[0]+twin_delta_pos[0], twin_sphere.center[1], twin_sphere.center[2]};
            twin_spheres.push_back(twin_sphere);
        }
        twin_ax2->set_spheres(twin_spheres);
        twin_axons.push_back(twin_ax2);
    }
}

bool AxonGammaDistribution::isSphereColliding(Dynamic_Sphere sph){
    
    Vector3d position = sph.center;
    double distance_to_be_inside = sph.max_radius + 10*barrier_tickness;

    std::vector<int> col_sphere_ids;

    for (unsigned i = 0; i < axons.size() ; i++){

        bool isinside = axons[i].isPosInsideAxon(position, distance_to_be_inside, true, col_sphere_ids);
        
        if (isinside){
            return true;
            break;
        }
    }
    return false;
}

bool AxonGammaDistribution::isAxonColliding(Axon* ax){
    for (unsigned i = 0; i < ax->spheres.size() ; i++){
        if (isSphereColliding(ax->spheres[i])){
            return true;
            break;
        } 
    }  
    return false;
}  

bool AxonGammaDistribution::find_next_center(Axon* ax, Dynamic_Sphere& s, vector<Eigen::Vector3d> centers, double dist_, double& rad, ostream& out, bool& collides_with_border, Eigen::Vector2d& twin_delta_pos,  std::vector<Eigen::Vector2d> phi_theta_colliding, int max_tries = 1000){

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

        tie(phi_to_target, theta_to_target) = phi_theta_to_target (curr_pos, ax->end,  out);

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
        if (acos(v1.dot(v2)) > M_PI/2){
            tries += 1;
            continue;
        }
        
        

        Dynamic_Sphere sphere_ (new_pos, rad,ax->volume_inc_perc,ax->swell, ax->id, 1, ax->active_state);
        
        // if sphere doesn't collide with environment
        if (!isSphereColliding(sphere_)){
            // if sphere collides with border planes
            if (check_borders(sphere_.center, sphere_.max_radius, twin_delta_pos)){
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

bool AxonGammaDistribution::fiber_collapse(std::vector<Eigen::Vector3d>& centers, int fibre_collapsed_nbr, ostream& out){
    
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



void AxonGammaDistribution::fill_with_spheres(Axon* ax, std::vector<Dynamic_Sphere>& spheres_to_add, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii){

    for (unsigned i = 1; i < centers.size(); i++){

        Dynamic_Sphere sphere (centers[i-1], sph_radii[i-1] ,ax->volume_inc_perc, ax->swell, ax->id, 1, ax->active_state, (i-1)*4);  
        spheres_to_add.push_back(sphere);
        // add spheres in between 
        for (unsigned j = 0; j < 3; j++){
            // center
            Eigen::Vector3d center_between = {(3-j)*centers[i-1][0]/4+(j+1)*centers[i][0]/4, (3-j)*centers[i-1][1]/4+(j+1)*centers[i][1]/4, (3-j)*centers[i-1][2]/4+(j+1)*centers[i][2]/4};
            // radius
            double radius_in_between = (3-j)*sph_radii[i-1]/4+(j+1)*sph_radii[i]/4;
            // sphere
            Dynamic_Sphere sphere_in_between (center_between, radius_in_between,ax->volume_inc_perc, ax->swell, ax->id, 1, ax->active_state, (i-1)*4+j+1);  
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
                sphere_ = Dynamic_Sphere (center_between, rad,ax->volume_inc_perc, ax->swell, ax->id, 1, ax->active_state, sphere_.id);  
                shrink_perc += 0.02;
            }
            sphere_in_between = sphere_;

            bool last_sphere_add =false;
            if (!last_sphere_add){  
                if (sphere_in_between.center[2]<= max_limits[2]){
                    spheres_to_add.push_back(sphere_in_between);
    
                    if (i == centers.size()-1 && j== 2 ){
                        Dynamic_Sphere sphere (centers[i], sph_radii[i] ,ax->volume_inc_perc, ax->swell, ax->id, 1, ax->active_state, (i-1)*4+j+2);  
                        spheres_to_add.push_back(sphere);

                    } 
                }
                else{
                        // distances to border 
                        double dist_i_1= abs((centers[i-1] - max_limits).norm());
                        double dist_i = abs((centers[i] - max_limits).norm());
                        // vector from centers[i-1] to centers[i] 
                        Eigen::Vector3d vector = (centers[i] - centers[i-1]).normalized();
                        // centers[i-1][2] + lambda*vector[2] = max_limits[2]
                        double lambda_ = (max_limits[2]-centers[i-1][2])/vector[2];
                        Eigen::Vector3d last_center = centers[i-1] +lambda_*vector;
                        double weight_i_1 = abs((last_center-centers[i]).norm())/abs((centers[i-1]-centers[i]).norm());
                        double weight_i = abs((last_center-centers[i-1]).norm())/abs((centers[i-1]-centers[i]).norm());
                        double last_radius = weight_i*sph_radii[i]+weight_i_1*sph_radii[i-1];
                        Dynamic_Sphere last_sphere(last_center, last_radius,ax->volume_inc_perc,ax->swell, ax->id, 1, ax->active_state, (i-1)*4+j+1);
                        last_sphere_add = true;
                        if (!isSphereColliding(last_sphere)){
                            spheres_to_add.push_back(last_sphere);
                        }
                }
            } 
        } 
    } 
} 


void AxonGammaDistribution::add_next_sphere(Dynamic_Sphere added_sphere, std::vector<Eigen::Vector3d>& centers, std::vector<double>& sph_radii){
    
    centers.push_back(added_sphere.center);
    sph_radii.push_back(added_sphere.min_radius);
        
}

void AxonGammaDistribution::shrink_sphere_rad(double& rad, double axon_rad, double& shrink_perc){
    rad = axon_rad *(1-shrink_perc);

}

void AxonGammaDistribution::find_shrinking_dichotomy(double& rad, double axon_rad, double min_shrink, double max_shrink, double& half_shrink){

    half_shrink = (max_shrink-min_shrink)/2+min_shrink;
    rad = axon_rad *(1-half_shrink);
    if (rad < 0.05*1e-3) {
        rad = 0.05*1e-3;
    }

}


std::vector<Dynamic_Sphere> AxonGammaDistribution::GrowAxon(Axon* ax, double distance_to_be_inside, bool& axon_collides_with_border, std::vector<Eigen::Vector2d>& all_twin_delta_pos, bool can_shrink, ostream& out){
    
    // list of all positions of the spheres add to axon in the growth
    std::vector<Eigen::Vector3d> centers;
    // radii of all spheres add to the centers list
    std::vector<double> sph_radii;
    // list of spheres to add
    std::vector<Dynamic_Sphere> spheres_to_add;

    // maximum of shrinking percentage
    double max_shrinking = 0.8;
    double min_radius_shrinking = ax->min_radius*(1-max_shrinking);
    // minimum radius of axon
    if (min_radius_shrinking < 0.05*1e-3) {
        min_radius_shrinking = 0.05*1e-3;
    }
    // managed to add 4 spheres
    bool spheres_added = false;
    // fiber collapsing is possible
    bool can_fiber_collapse = true;
    int fibre_collapsed_nbr = 0;
    // if axon collides with border plane
    axon_collides_with_border = false;
    // phi and thetas to exclude in tries
    std::vector<std::vector<Eigen::Vector2d>> all_phi_theta_colliding;

    // first sphere to add
    Dynamic_Sphere s1(ax->begin, ax->min_radius,ax->volume_inc_perc,ax->swell, ax->id, 1, ax->active_state);

    if(isSphereColliding(s1)) {
        return spheres_to_add;
    }
    else{
        // if the first sphere touches boundaries
        Eigen::Vector2d twin_delta_pos_;
        if(check_borders(s1.center, s1.max_radius, twin_delta_pos_)){
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
                return spheres_to_add;
            }

        }
        else{
            add_next_sphere(s1, centers, sph_radii);
        }
    }
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
        spheres_added = find_next_center(ax, sphere_to_add, centers, dist_, rad, out, collides_with_border, twin_delta_pos, phi_theta_colliding);
        
        if (spheres_added){
            add_next_sphere(sphere_to_add, centers, sph_radii);
            // phi and theta that didn't work are saved 
            all_phi_theta_colliding.push_back(phi_theta_colliding);
            if(collides_with_border){
                axon_collides_with_border = true;
                // if all_twin_delta_pos is empty or twin_delta_pos not already in all_twin_delta_pos
                if (all_twin_delta_pos.size() == 0 || std::find(all_twin_delta_pos.begin(), all_twin_delta_pos.end(),twin_delta_pos) ==all_twin_delta_pos.end()) {
                    all_twin_delta_pos.push_back(twin_delta_pos);
                }
            }
            out << "axon_collides_with_border "<< axon_collides_with_border << endl;
            //out << "managed to add spheres :" << spheres_added<< endl;
            out << "centers size :" << centers.size() << endl;
            out << "all_phi_theta_colliding size :" << all_phi_theta_colliding.size() << endl;
        }

        else{

            // fiber collapse 
            // check and can do fiber collapsing
            can_fiber_collapse = fiber_collapse(centers, fibre_collapsed_nbr, out);
            if (can_fiber_collapse && all_phi_theta_colliding.size() > 0){
                std::vector<Eigen::Vector2d> phi_theta_colliding_collapse;
                out << "fiber collapse"<< endl;
                // delete all necessary phi and thetas from previous steps
                all_phi_theta_colliding.pop_back();
                // find next center while not retaking phis and thetas that didn't work
                phi_theta_colliding_collapse = all_phi_theta_colliding[all_phi_theta_colliding.size()-1]; 
                // add next center
                spheres_added = find_next_center(ax, sphere_to_add, centers, dist_, rad, out, collides_with_border, twin_delta_pos, phi_theta_colliding_collapse);
                if (spheres_added){
                    add_next_sphere(sphere_to_add, centers, sph_radii);
                    // add thetas and phis that didn't work in list
                    all_phi_theta_colliding.push_back(phi_theta_colliding_collapse);
                    // border
                    if(collides_with_border){
                        axon_collides_with_border = true;
                        if (all_twin_delta_pos.size() == 0 || std::find(all_twin_delta_pos.begin(), all_twin_delta_pos.end(),twin_delta_pos) ==all_twin_delta_pos.end()) {
                            all_twin_delta_pos.push_back(twin_delta_pos);
                        }
                    }
                }
                fibre_collapsed_nbr += 1;
                out << "managed to add spheres after fiber collapse :" << spheres_added << endl;
                out << "centers size :" << centers.size() << endl;
            }
            // if cannot fiber collapse -> shrinking 
            else if (can_shrink){ 
                std::vector<Eigen::Vector2d> phi_theta_colliding_shrink;
                Dynamic_Sphere sphere_to_add_;
                int max_tries = 10000;
                phi_theta_colliding_shrink.clear();
                bool can_add_with_max_shrinking = find_next_center(ax, sphere_to_add_, centers, dist_, min_radius_shrinking, out, collides_with_border, twin_delta_pos, phi_theta_colliding_shrink, max_tries);
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
                                axon_collides_with_border = true;
                                // if all_twin_delta_pos is empty or twin_delta_pos not already in all_twin_delta_pos
                                if (all_twin_delta_pos.size() == 0 || std::find(all_twin_delta_pos.begin(), all_twin_delta_pos.end(),twin_delta_pos) ==all_twin_delta_pos.end()) {
                                    all_twin_delta_pos.push_back(twin_delta_pos);
                                }
                            }
                            achieved = true;
                            out << "    sphere shrinks to :"<< rad << endl;
                            out << "centers size :" << centers.size() << endl;
                            out << "all_phi_theta_colliding size :" << all_phi_theta_colliding.size() << endl;
                            }
                        else{   
                            // rad = shrinked radius
                            find_shrinking_dichotomy(rad, ax->min_radius, lower_boundary, upper_boundary, middle_boundary);
                            // distance between spheres is the maximum between the radius of the two
                            dist_ = max(rad, sph_radii[sph_radii.size()-1]);
                            // try shrinking at half 
                            phi_theta_colliding_shrink.clear();
                            dichotomy_check = find_next_center(ax, sphere_to_add_, centers, dist_, rad, out, collides_with_border, twin_delta_pos, phi_theta_colliding_shrink);
                                
                            out << "    checking shrink : " << middle_boundary  << "%, result :"<< dichotomy_check << endl;
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
                    return spheres_to_add;
                    break;

                }  

            }
            else{
                // this will be empty
                out << "cannot add max shrinked sphere"<< endl;
                return spheres_to_add;
                break;
            }
                
        }


    out << "last center : (" << centers[centers.size()-1][0] << ","<< centers[centers.size()-1][1] << ","<< centers[centers.size()-1][2] <<"), ax->end[2] :" << ax->end[2] << endl;
    //out << "small voxel size :" << small_voxel_size << endl;
    }while (centers[centers.size()-1][2] < ax->end[2]);
    
    fill_with_spheres(ax, spheres_to_add,centers, sph_radii);



    return spheres_to_add;

}   