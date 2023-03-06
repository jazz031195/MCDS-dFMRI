#include "obstacle.h"
#include <math.h>

Obstacle::Obstacle():id(-1),percolation(0),T2(0)
{
}

bool Obstacle::checkCollision(Walker& walker, Eigen::Vector3d const& step, double const& step_lenght, Collision& colision)
{
    return false;
}

void Obstacle::elasticBounceAgainsPlane(Eigen::Vector3d const& ray_origin, Eigen::Vector3d const& normal, double const& step_length, Eigen::Vector3d &step_dir) const
{

    Eigen::Vector3d ray =  -step_dir;//
    double rn = ray.dot(normal);

    // Caso 3) ni cerca ni paralela
    step_dir = -ray + 2.0*normal*rn;

    //step = (rn>0.0)?normal:(-normal);

}

double Obstacle::minDistance(Walker const& walker) const
{
    return 0;

}
