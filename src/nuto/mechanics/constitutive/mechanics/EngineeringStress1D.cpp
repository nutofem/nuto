// $Id$
#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"

// constructor
NuTo::EngineeringStress1D::EngineeringStress1D()
{
    this->mEngineeringStress = 0.0;
}

// number of components
unsigned int NuTo::EngineeringStress1D::GetNumberOfComponents() const
{
    return 1;
}

// get Engineering stress
const double* NuTo::EngineeringStress1D::GetData() const
{
    return &mEngineeringStress;
}

// info routine
void NuTo::EngineeringStress1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): " << this->mEngineeringStress << std::endl;
}
