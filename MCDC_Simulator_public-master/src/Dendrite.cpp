#include "Dendrite.h"

#include "Eigen/Dense"
#include <iostream>

#include <cmath>
#include <typeinfo>

using namespace Eigen;
using namespace std;

int Dendrite::nb_dendrites = 0;

Dendrite::Dendrite()
{
    id          = nb_dendrites++;
    dendrite_id = -1;
}

Dendrite::Dendrite(std::vector<Axon> const& subbranches_):Dendrite()
{
    subbranches = subbranches_;
}

Dendrite::Dendrite(Dendrite const& dendrite)
{
    subbranches = dendrite.subbranches;
}

bool Dendrite::checkCollision(Walker &walker, Vector3d const&step_dir, double const&step_lenght, Collision &collision)
{
    // TO IMPLEMENT
    return false;
}

bool Dendrite::isPosInsideDendrite(Vector3d const& position,  double const& barrier_thickness, bool const& swell_, int& in_soma_index, int& in_dendrite_index)
{
    // TO IMPLEMENT
    return false;
}

double Dendrite::minDistance(Vector3d const& pos) const
{
    // TO IMPLEMENT
    return 0;
}

double Dendrite::minDistance(Walker const& walker) const
{
    // TO IMPLEMENT
    return 0;
}

void Dendrite::add_subbranch(Axon const& subbranch)
{
    subbranches.push_back(subbranch);
}

void Dendrite::set_dendrite(vector<Axon> const& subbranches_)
{
    subbranches = subbranches_;
}

int Dendrite::get_nb_subbranches()
{
    return subbranches.size();
}

double Dendrite::volumeDendrite() const
{
    double volume = 0;
    for (size_t i=0; i < subbranches.size(); i++)
        volume += subbranches[i].volumeAxon();
    
    return volume;
}