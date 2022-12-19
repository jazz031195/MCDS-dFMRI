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

void AxonGammaDistribution::createGammaSubstrate()
{
    // generate the gamma distribution
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    uint repetition = 40;
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
    std::vector<bool> bool_swell_ax_id(number_swelling_axons, false);

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

    uint adjustments = 0;
    // We increease 1% the total area. (Is prefered to fit all the cylinders than achieve a perfect ICVF.)
    double adj_increase = icvf * 0.01;


    while (!achieved)
    {

        double target_icvf = this->icvf + adjustments * adj_increase;
        // computes max_limits
        computeMinimalSize(radiis, target_icvf, max_limits);


        for (uint t = 0; t < repetition; t++)
        {
            vector<Axon> axons_to_add;

            axons.clear();
            for (unsigned i = 0; i < num_obstacles; i++)
            {
                
                unsigned stuck = 0;

                while (++stuck <= 1000)
                {

                    double t = udist(gen);
                    double x = (t * max_limits[0]) + (1 - t) * min_limits[0];
                    t = udist(gen);
                    double y = (t * max_limits[1] + (1 - t) * min_limits[1]);
                    double z = 0;

                    Vector3d Q = {x, y, z + min_limits[2]};
                    Vector3d D = {x, y, z + max_limits[2]};


                    Axon ax(radiis[i], Q, D, volume_inc_perc, bool_swell_ax_id[i], active_state);
                

                    double min_distance;
                    // creates dyn_cylinders_to_add
                    bool collision = checkForCollition(ax, min_limits, max_limits, axons_to_add, min_distance);
                    
                    if (!collision)
                    {
                        for (unsigned j = 0; j < axons_to_add.size(); j++){
                            axons.push_back(axons_to_add[j]);
                        }
                        break;
                    }
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
            break;
        }
    }

    axons = best_axons;
    max_limits = best_max_limits;

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
        out << " " << axons[i].begin[0] * 1e3  << " " << axons[i].begin[1] * 1e3  << " " 
        << axons[i].begin[2] * 1e3  << " " << axons[i].end[0] * 1e3  << " " 
        << axons[i].end[1] * 1e3  <<" "<< axons[i].end[2] * 1e3  <<" " 
        << axons[i].radius* 1e3 << " " << axons[i].swell 
        << endl;

    }
}

bool AxonGammaDistribution::checkForCollition(Axon ax, Vector3d min_limits, Vector3d max_limits, std::vector<Axon>& axons_to_add, double &min_distance)
{

    axons_to_add.clear();
    checkBoundaryConditions(ax, axons_to_add, min_limits, max_limits);

    min_distance = 1e10;

    bool collision = false;


    for (unsigned i = 0; i < axons.size(); i++)
    {
        for (unsigned j = 0; j < axons_to_add.size(); j++)
        {

            double distance = (axons[i].begin - axons_to_add[j].begin).norm();
            // give enough space so that each cylinder can potentially swell
            if (distance - (axons[i].radius + axons_to_add[j].radius) < 1e-15)
            {
                min_distance = 0;
                collision = true;
                break;
            }
            if (distance < min_distance)
                min_distance = distance;
        }
    }


    // we need to check that the cylinders to add don't interse4ct each other (very small voxel sizes)
    for (unsigned i = 0; i < axons_to_add.size() - 1; i++)
    {

        for (unsigned j = i + 1; j < axons_to_add.size(); j++)
        {

            double distance = (axons_to_add[i].begin  - axons_to_add[j].begin ).norm();
            // give enough space so that each cylinder can potentially swell

            if (distance - (axons_to_add[i].radius + axons_to_add[j].radius) < 1e-15)
            {
                min_distance = 0;
                collision = true;
                break;
            }
            else if (distance < min_distance)
                min_distance = distance;
            
        }
    }


    return collision;
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

    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]);

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
        AreaC += M_PI * rad * rad;
        num_no_repeat++;
    }
    return AreaC / AreaV;
}

void AxonGammaDistribution::checkBoundaryConditions(Axon ax, std::vector<Axon> &axons_to_add, Vector3d min_limits, Vector3d max_limits)
{
    vector<Axon> to_add;
    
    to_add.push_back(ax);


    for (int i = 0; i < 2; i++)
    {

        double rad = ax.radius;

        Vector3d P = ax.begin;
        Vector3d Q = ax.end;


        if (ax.begin[i] + rad >= max_limits[i])
        {
            
            P[i]  = ax.begin[i] + min_limits[i] - max_limits[i];
            Q[i]=  ax.end[i] + min_limits[i] - max_limits[i];
            Axon tmp (ax.radius, P, Q, ax.volume_inc_perc, ax.active_state);
            to_add.push_back(tmp);

        }


        if (ax.begin[i] - rad <= min_limits[i])
        {
            P[i]  = ax.begin[i] - min_limits[i] + max_limits[i];
            Q[i]=  ax.end[i] - min_limits[i] + max_limits[i];
            Axon tmp (ax.radius, P, Q, ax.volume_inc_perc, ax.active_state);
            to_add.push_back(tmp);

        }
    }



    if (to_add.size() == 3)

        for (unsigned j = 1; j < 3; j++)
        {
            Axon jkr(to_add[j]);


            for (int i = 0; i < 2; i++)
            {
                double rad = jkr.radius;

                Vector3d P = jkr.begin;
                Vector3d Q = jkr.end;

                if (jkr.begin[i] + rad >= max_limits[i])
                {

                    P[i]  = jkr.begin[i] + min_limits[i] - max_limits[i];
                    Q[i]=  jkr.end[i] + min_limits[i] - max_limits[i];
                    Axon tmp (jkr.radius, P, Q, jkr.volume_inc_perc, jkr.active_state);
                    to_add.push_back(tmp);


                }

                if (jkr.begin[i] - rad <= min_limits[i])
                {

                    P[i]  = jkr.begin[i] - min_limits[i] + max_limits[i];
                    Q[i]=  jkr.end[i] - min_limits[i] + max_limits[i];
                    Axon tmp (jkr.radius, P, Q, jkr.volume_inc_perc, jkr.active_state);
                    to_add.push_back(tmp);


                }
            }
        }

    for (unsigned i = 0; i < to_add.size(); i++)
    {
        bool rep = false;
        for (unsigned j = 0; j < axons_to_add.size(); j++)
        {
            if ((fabs(to_add[i].begin[0] - axons_to_add[j].begin[0]) < 1e-12) && (fabs(to_add[i].begin[1] - axons_to_add[j].begin[1]) < 1e-12))
            {
                rep = true;
                break;
            }
        }
        

        if (rep == false)
        {
            axons_to_add.push_back(to_add[i]);

        }

    }

}
