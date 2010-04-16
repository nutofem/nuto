// $Id$
#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"

// constructor
NuTo::EngineeringStress2D::EngineeringStress2D()
{
    this->mEngineeringStress[0] = 0.0;
    this->mEngineeringStress[1] = 0.0;
    this->mEngineeringStress[2] = 0.0;
}

// number of components
unsigned int NuTo::EngineeringStress2D::GetNumberOfComponents() const
{
    return 3;
}

// get Engineering stress
const double* NuTo::EngineeringStress2D::GetData() const
{
    return this->mEngineeringStress;
}

// info routine
void NuTo::EngineeringStress2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): " << this->mEngineeringStress[0] << ", "
              << this->mEngineeringStress[1] << ", " << this->mEngineeringStress[2]  << std::endl;
}
