#include "dyncylindergammadistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;

DynCylinderGammaDistribution::DynCylinderGammaDistribution(double dyn_perc_, double volume_inc_perc_, unsigned num_cyl, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l, float min_radius_)
{
    dyn_perc = dyn_perc_;
    volume_inc_perc = volume_inc_perc_;
    num_obstacles = num_cyl;
    alpha = a;
    beta  = b;
    icvf = icvf_;
    min_limits = min_l;
    max_limits = max_l;
    min_radius = min_radius_;
    icvf_current = 0;
    dyn_cylinders.clear();

}
void DynCylinderGammaDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d &l)
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

void DynCylinderGammaDistribution::displayGammaDistribution()
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
        SimErrno::info(message,std::cout);
    }
    message = "9-10:" + std::string(p[9]*nstars/nrolls,'*') ;
    SimErrno::info(message,std::cout);
    message = ">10: " +  std::string(p[10]*nstars/nrolls,'*') + "\n" ;
    SimErrno::info(message,std::cout);
}

void DynCylinderGammaDistribution::createGammaSubstrate()
{
    // generate the gamma distribution
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha, beta);
    std::vector<Dynamic_Cylinder*> dyn_cylinders_to_add;
    Eigen::Vector3d best_max_limits;
    min_limits = {0., 0., 0.};

    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0, 1);
    std::vector<double> radiis(num_obstacles, 0);

    int number_swelling_cylinders = int(num_obstacles * dyn_perc);
    std::vector<bool> bool_swell_cyl_id(num_obstacles, false);
    int perc_;

    int tried = 0;


    // create list of ids that correspond to dynamic cylinders based on percentage input
    for (unsigned i = 0; i < number_swelling_cylinders; ++i)
    {
        int random_id = rand() % num_obstacles;
        if (i > 0)
        {
            while (bool_swell_cyl_id[random_id]) {
                random_id = rand() % num_obstacles;
            }
        }
        bool_swell_cyl_id[random_id] = true;
    }

    for (unsigned i = 0; i < num_obstacles; ++i)
    {

        if (tried > 10000)
        {
            string message = " Radii distribution cannot be sampled [Min. radius Error]\n";
            SimErrno::error(message, std::cout);
            assert(0);
        }
        double jkr = distribution(generator);

        if (jkr < min_radius  )
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

    computeMinimalSize(radiis, icvf, max_limits);

    for (unsigned i = 0; i < num_obstacles; i++)
    {
                
        unsigned stuck = 0;

        bool cyl_added = false;
        while (++stuck <= 1000 && !cyl_added)
        {


            double t = udist(gen);
            double x = (t * max_limits[0]) + (1 - t) * min_limits[0];
            t = udist(gen);
            double y = (t * max_limits[1] + (1 - t) * min_limits[1]);
            double z = 0;

            Vector3d Q = {x, y, z};
            Vector3d D = {x, y, z + 1};

            Dynamic_Cylinder* cyl = new Dynamic_Cylinder ( dyn_cylinders.size(), Q, D, volume_inc_perc, bool_swell_cyl_id[i], radiis[i]);

            double min_distance;
            // creates dyn_cylinders_to_add

            bool collision = checkForCollition(cyl, min_limits, max_limits, dyn_cylinders_to_add, min_distance);

            if (!collision)
            {
                for (unsigned j = 0; j < dyn_cylinders_to_add.size(); j++){
                    Dynamic_Cylinder* cyl_to_add = dyn_cylinders_to_add[j];
                    cyl_to_add->id = dyn_cylinders.size();
                    dyn_cylinders.push_back(cyl_to_add);
                }
                cyl_added= true;

            }
                    
        }

    icvf_current = computeICVF(dyn_cylinders, min_limits, max_limits, perc_);
    std::cout << "icvf_current :" << icvf_current << endl;


    } // end for cylinders

    string message = "Percentage of cylinders selected: " + to_string(double(perc_) / radiis.size() * 100.0) + "%,\nICVF achieved: " + to_string(icvf_current * 100) + "  (" + to_string(int((icvf_current / icvf * 100))) + "% of the desired icvf)\n";
    SimErrno::info(message, std::cout);
}

void DynCylinderGammaDistribution::printSubstrate(ostream &out)
{
    double scale = 1;
    out << scale << endl;
    out << volume_inc_perc << endl;
    out << dyn_perc << endl;
    out << icvf_current << endl;
    out << min_limits[2] << endl;
    out << max_limits[2] << endl;

    for (unsigned i = 0; i < dyn_cylinders.size(); i++)
    {
        
        out <<  dyn_cylinders[i]->P[0]/scale << " " << dyn_cylinders[i]->P[1]/scale << " " << dyn_cylinders[i]->P[2]/scale << " "
            << dyn_cylinders[i]->Q[0]/scale << " " << dyn_cylinders[i]->Q[1]/scale<< " " << dyn_cylinders[i]->Q[2]/scale << " "
            << dyn_cylinders[i]->min_radius/scale  << " " << dyn_cylinders[i]->swell << endl;
    }
}

bool DynCylinderGammaDistribution::checkForCollition(Dynamic_Cylinder* cyl, Vector3d min_limits, Vector3d max_limits, std::vector<Dynamic_Cylinder*>& dyn_cylinders_to_add, double &min_distance)
{

    dyn_cylinders_to_add.clear();
    checkBoundaryConditions(cyl, dyn_cylinders_to_add, min_limits, max_limits);

    double security_distance = 2*barrier_tickness;

    min_distance = 1e10;

    bool collision = false;

    for (unsigned i = 0; i < dyn_cylinders.size(); i++)
    {
        for (unsigned j = 0; j < dyn_cylinders_to_add.size(); j++)
        {
            double distance = (dyn_cylinders[i]->P - dyn_cylinders_to_add[j]->P).norm();
            // give enough space so that each cylinder can potentially swell
            if (distance - (dyn_cylinders[i]->radius + dyn_cylinders_to_add[j]->radius) < security_distance)
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

    for (unsigned i = 0; i < dyn_cylinders_to_add.size() - 1; i++)
    {
        if (dyn_cylinders_to_add.size() > 0) {
            for (unsigned j = i + 1; j < dyn_cylinders_to_add.size(); j++)
            {

                double distance = (dyn_cylinders_to_add[i]->P - dyn_cylinders_to_add[j]->P).norm();
                // give enough space so that each cylinder can potentially swell
                if (distance - (dyn_cylinders_to_add[i]->radius + dyn_cylinders_to_add[j]->radius) < 1e-15)
                {
                    min_distance = 0;
                    collision = true;
                    break;
                }
                else if (distance < min_distance)
                    min_distance = distance;
                
            }
        }
    }


    return collision;
}

/*
WARNING: The way we discard repeated cylinders is using radius. Repreated radius (like really the same)
are considered the same. This becasuse we don't track wich cylinders had to be replciated to mantain the voxel
symmetry
*/

double DynCylinderGammaDistribution::computeICVF(std::vector<Dynamic_Cylinder*> &dyn_cylinders, Vector3d &min_limits, Vector3d &max_limits, int &num_no_repeat)
{

    if (dyn_cylinders.size() == 0)
        return 0;

    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]);

    double AreaC = 0;

    // using a lambda function:
    //std::sort(dyn_cylinders.begin(), dyn_cylinders.end(), [](const Dynamic_Cylinder a, Dynamic_Cylinder b) -> bool
    //          { return a.radius > b.radius; });

    double rad_holder = -1;
    num_no_repeat = 0;
    for (uint i = 0; i < dyn_cylinders.size(); i++)
    {

        if (fabs(rad_holder - dyn_cylinders[i]->radius) < 1e-15)
        {
            continue;
        }
        else
        {
            rad_holder = dyn_cylinders[i]->radius;
        }

        double rad = dyn_cylinders[i]->radius;
        AreaC += M_PI * rad * rad;
        num_no_repeat++;
    }
    return AreaC / AreaV;
}

void DynCylinderGammaDistribution::checkBoundaryConditions(Dynamic_Cylinder* cyl, std::vector<Dynamic_Cylinder*> &dyn_cylinders_to_add, Vector3d min_limits, Vector3d max_limits)
{
    vector<Dynamic_Cylinder*> to_add;
    
    to_add.push_back(cyl);


    for (int i = 0; i < 2; i++)
    {

        double rad = cyl->radius;

        if (cyl->P[i] + rad >= max_limits[i])
        {

            Dynamic_Cylinder* tmp = cyl;
            tmp->P[i] += min_limits[i] - max_limits[i];
            tmp->Q[i] += min_limits[i] - max_limits[i];
            to_add.push_back(tmp);


        }

        if (cyl->P[i] - rad <= min_limits[i])
        {
            Dynamic_Cylinder* tmp = cyl;
            tmp->P[i] += max_limits[i] - min_limits[i];
            tmp->Q[i] += max_limits[i] - min_limits[i];
            to_add.push_back(tmp);
        }
    }

    if (to_add.size() == 3)
        for (unsigned j = 1; j < 3; j++)
        {
            Dynamic_Cylinder* jkr = to_add[j];
            for (int i = 0; i < 2; i++)
            {
                double rad = cyl->radius;

                if (jkr->P[i] + rad >= max_limits[i])
                {
  
                    Dynamic_Cylinder* tmp = jkr;
                    tmp->P[i] += min_limits[i] - max_limits[i];
                    tmp->Q[i] += min_limits[i] - max_limits[i];
                    to_add.push_back(tmp);


                }

                if (jkr->P[i] - rad <= min_limits[i])
                {
                    Dynamic_Cylinder* tmp = jkr;
                    tmp->P[i] += max_limits[i] - min_limits[i];
                    tmp->Q[i] += max_limits[i] - min_limits[i];
                    to_add.push_back(tmp);

                }
            }
  
        }

    for (unsigned i = 0; i < to_add.size(); i++)
    {
        bool rep = false;
        for (unsigned j = 0; j < dyn_cylinders_to_add.size(); j++)
        {
            if ((fabs(to_add[i]->P[0] - dyn_cylinders_to_add[j]->P[0]) < 1e-12) && (fabs(to_add[i]->P[1] - dyn_cylinders_to_add[j]->P[1]) < 1e-12))
            {
                rep = true;
                break;
            }
        }
        

        if (rep == false)
        {
            dyn_cylinders_to_add.push_back(to_add[i]);

        }

    }

}
