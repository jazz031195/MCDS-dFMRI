//!  Cylinder Obstacle Derived Class =============================================================/
/*!
*   \details   Cylinder class derived from an Obstacle. Defines infinite long cylinders
*              in the direction set by P,Q.
*   \author    Jasmine Nguyen-Duc
*   \date      September 2022 
*   \version   1.42
=================================================================================================*/


#ifndef AXON_H
#define AXON_H

#include "dynamic_sphere.h"


/// @brief 
class Axon : public Obstacle
{
public:

    static int count;
    std::vector<Dynamic_Sphere> spheres; 
    std::vector<double> radii;
    std::vector<Eigen::Vector3d> centers; 
    bool swell;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();


    Axon(std::vector<double> radii_, std::vector<Eigen::Vector3d> centers_, bool swell_ = false):radii(radii_),centers(centers_),swell(swell_){
        
        for (unsigned i=0; i< radii.size(); ++i){
            Dynamic_Sphere sphere(centers_[i], radii_[i]);  
            spheres[i] = sphere;
        }
        id = count++;
    }
    Axon(Axon const &ax);



};

#endif // AXON_H
