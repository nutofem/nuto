// $Id$
#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"

// constructor
NuTo::EngineeringStress3D::EngineeringStress3D()
{
    for (unsigned int count = 0; count < 6; count++)
    {
        this->mEngineeringStress[count] = 0.0;
    }
}

// number of components
unsigned int NuTo::EngineeringStress3D::GetNumberOfComponents() const
{
    return 6;
}

// get Engineering stress
const double* NuTo::EngineeringStress3D::GetData() const
{
    return this->mEngineeringStress;
}


// info routine
void NuTo::EngineeringStress3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): "
              << this->mEngineeringStress[0] << ", " << this->mEngineeringStress[1] << ", " << this->mEngineeringStress[2] << ", "
              << this->mEngineeringStress[3] << ", " << this->mEngineeringStress[4] << ", " << this->mEngineeringStress[5] << std::endl;
}
