// $Id$

#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress3D.h"

// constructor
NuTo::SecondPiolaKirchhoffStress3D::SecondPiolaKirchhoffStress3D()
{
    for (unsigned int count = 0; count < 6; count++)
    {
        this->mSecondPiolaKirchhoffStress[count] = 0.0;
    }
}

// number of components
unsigned int NuTo::SecondPiolaKirchhoffStress3D::GetNumberOfComponents() const
{
    return 6;
}

// get second Piola-Kirchhoff stress
const double* NuTo::SecondPiolaKirchhoffStress3D::GetData() const
{
    return this->mSecondPiolaKirchhoffStress;
}

// info routine
void NuTo::SecondPiolaKirchhoffStress3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of second Piola-Kirchhoff stress tensor (vector notation): "
              << this->mSecondPiolaKirchhoffStress[0] << ", "
              << this->mSecondPiolaKirchhoffStress[1] << ", "
              << this->mSecondPiolaKirchhoffStress[2] << ", "
              << this->mSecondPiolaKirchhoffStress[3] << ", "
              << this->mSecondPiolaKirchhoffStress[4] << ", "
              << this->mSecondPiolaKirchhoffStress[5] << std::endl;
}

