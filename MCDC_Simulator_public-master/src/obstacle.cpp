#include "obstacle.h"
#include <math.h>

Obstacle::Obstacle():id(-1),percolation(0),T2(0)
{
}

bool Obstacle::checkCollision(Walker& walker, Eigen::Vector3d const& step, double const& step_lenght, Collision& colision)
{
    return false;
}

void Obstacle::elasticBounceAgainsPlane(Eigen::Vector3d& ray_origin, Eigen::Vector3d& normal, double& step_length, Eigen::Vector3d &step_dir) 
{

    Eigen::Vector3d ray =  -step_dir;//
    double rn = ray.dot(normal);

    // Caso 3) ni cerca ni paralela
    step_dir = -ray + 2.0*normal*rn;

    //step = (rn>0.0)?normal:(-normal);

}

void Obstacle::elasticBounceAgainsPlane_intra(Eigen::Vector3d &ray_origin, Eigen::Vector3d &normal, double &t, Eigen::Vector3d &step)
{

    Eigen::Vector3d ray =  (-t*step).normalized();//

    double rn = ray.dot(normal);
    if (cos(rn)== 0){
        step = normal;
    } 
    else{
        step = -ray + 2.0*normal*rn;
    } 


}

void Obstacle::elasticBounceAgainsPlane_extra(Eigen::Vector3d &ray_origin, Eigen::Vector3d &normal, double &t, Eigen::Vector3d &step)
{

    Eigen::Vector3d ray =  (-t*step).normalized();//

    double rn = ray.dot(normal);
    if (cos(rn) <= 0){
        step = -ray + 2.0*normal*rn;
    } 


}

double Obstacle::minDistance(Walker &w)
{
    return 0;
}
