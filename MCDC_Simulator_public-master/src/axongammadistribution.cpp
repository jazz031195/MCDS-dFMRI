#include "axongammadistribution.h"
#include "grow_axons.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"
#include <chrono>



using namespace std;
using namespace Eigen;
using namespace std::chrono;

//* Auxiliare method to split words in a line using the spaces*//
template<typename Out>
void _split_(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> _split_line(const std::string &s, char delim) {
    std::vector<std::string> elems;
    _split_(s, delim, std::back_inserter(elems));
    return elems;
}

AxonGammaDistribution::AxonGammaDistribution (double dyn_perc_,double volume_inc_perc_, unsigned num_ax, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l, double min_radius, bool active_state_, double c2_, bool tortuous_, double step_length_, string gamma_from_file_)
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
    gamma_from_file = gamma_from_file_;

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

/*
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
*/

/*
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
                Axon ax = Axon (axons[i]);
                ax.begin = new_begin;
                ax.end = new_end;
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
                    ax.set_spheres(new_spheres);
                    axons_to_add.push_back(ax);
                }

            }
        }
    }
    // extend axons with axons to add
    axons.reserve(axons.size() + distance(axons_to_add.begin(),axons_to_add.end()));
    axons.insert(axons.end(),axons_to_add.begin(),axons_to_add.end());
}
*/

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
    // read gamma distribution from file
    if ( beta== 0.0 && alpha == 0.0){
        std::ifstream in(gamma_from_file);
        if(!in){
            std::cout <<  "[ERROR] Unable to open:" << gamma_from_file << std::endl;
            return;
        }
        int i = 0;
        for( std::string line; getline( in, line ); )
        {
            if (i< num_obstacles) {
                double jkr = stod(_split_line(line,' ')[0]);
                if (jkr < 0.05){
                    jkr = 0.05;
                }
                radiis[i] = jkr*1e-3; // WE CONVERT FROM UM TO MM HERE
            }
            i += 1;
        }
    }
    else{
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
    }

    // sorts the radii 
    std::sort(radiis.begin(), radiis.end(), [](const double a, double b) -> bool
              { return a > b; });

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

        int number_axons_wo_twins = 0;


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
                int id_ = axons.size();
                Axon ax  = Axon (id_, radiis[i], Q, D, volume_inc_perc, active_state, bool_swell_ax_id[i], 1);

                Growth* growth1 = new Growth (ax, axons, small_voxel_size, tortuous);
                
                Growth* growth2 = new Growth (ax, axons, small_voxel_size, tortuous);
                
                std::vector<std::thread> growing_threads;
                
                // grow in parallel towards begin and end
                growing_threads.push_back(std::thread(&Growth::GrowInParallel,growth1, ax.end));
                growing_threads.push_back(std::thread(&Growth::GrowInParallel, growth2, ax.begin));

                for (unsigned int j =0; j < growing_threads.size(); j++){
                    growing_threads[j].join();
                }
                out << "growth1 axon :" << growth1->axon_to_grow.spheres.size() << endl;
                out << "growth2 axon :" << growth2->axon_to_grow.spheres.size()<< endl;

                // if didn't manage to grow both sides
                if(!(growth1->success) || !(growth2->success)){
                    out << "didn't manage" <<endl;

                    i -= 1; 
                    stuck += 1;
                    delete growth1;
                    delete growth2;
                    continue;
                }
                else{
                    combine_axons_and_save(growth1->axon_to_grow, growth2->axon_to_grow, axons.size());
                    number_axons_wo_twins += 2;
                    int min_nbr_twins;
                    if(growth1->twin_axons.size() > growth2->twin_axons.size()) {
                        min_nbr_twins = growth2->twin_axons.size();
                    }
                    else {
                        min_nbr_twins = growth1->twin_axons.size();
                    }
                    // add paired twins as one axon
                    for (unsigned j = 0; j < min_nbr_twins; j++){
                        // check if there are some pairs in twins to combine them
                        if (growth1->twin_axons[j].spheres[-1].center[2] == growth2->twin_axons[j].spheres[-1].center[2] ) {
                            combine_axons_and_save(growth1->twin_axons[j], growth2->twin_axons[j], axons.size());
                            // delete from list of twins
                            growth1->twin_axons.erase(growth1->twin_axons.begin() + j);
                            growth2->twin_axons.erase(growth2->twin_axons.begin() + j);
                        }
                    }
                    // add rest of the twins
                    for (unsigned j = 0; j < growth1->twin_axons.size(); j++){
                        growth1->twin_axons[j].id = axons.size();
                        axons.push_back(growth1->twin_axons[j]);
                    }
                    for (unsigned j = 0; j < growth2->twin_axons.size(); j++){
                        growth2->twin_axons[j].id = axons.size();
                        axons.push_back(growth2->twin_axons[j]);
                    }
                    delete growth1;
                    delete growth2;
                    out << "new axons size :" << axons.size()<<endl;
                } 
            }
            display_progress(number_axons_wo_twins, num_obstacles*2);
                    

        } // end for axons

        // check ICVF 
        icvf_current = computeICVF();
        achieved = true;
    }
    num_obstacles = number_axons_placed;

    std::sort(axons.begin(), axons.end(), [](const Axon a, Axon b) -> bool
              { return a.id < b.id; });


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

void AxonGammaDistribution::combine_axons_and_save(Axon axon1, Axon axon2, int id_) 
{
    Eigen::Vector3d D = axon1.spheres[axon1.spheres.size()-1].center;
    Eigen::Vector3d Q = axon2.spheres[axon2.spheres.size()-1].center;

    Axon* ax  = new Axon (id_, axon1.min_radius, D, Q, axon1.volume_inc_perc, axon1.active_state, axon1.swell);
    std::vector<Dynamic_Sphere> combined_spheres;
    combined_spheres = axon1.spheres;
    // spheres of axon2 without first element which is also in axon1
    std::vector<Dynamic_Sphere> combined_spheres2  (axon2.spheres.begin() + 1, axon2.spheres.end());
    // combine the spheres
    combined_spheres.reserve(combined_spheres.size() + distance(combined_spheres2 .begin(),combined_spheres2 .end()));
    combined_spheres.insert(combined_spheres.end(),combined_spheres2 .begin(),combined_spheres2 .end());

    // sort by z comp of center
    std::sort(combined_spheres.begin(),combined_spheres.end(), [](const Dynamic_Sphere a, Dynamic_Sphere b) -> bool
              { return a.center[2] < b.center[2]; });
              
    for (unsigned i = 0; i < combined_spheres.size(); i++) {
        combined_spheres[i].id = i;
        combined_spheres[i].ax_id = id_;
    }

    ax->set_spheres(combined_spheres);
    axons.push_back(*ax);

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

    double AreaC = 0;

    double tortuosity;

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
        else if (axons[i].spheres.size() > 1){

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
