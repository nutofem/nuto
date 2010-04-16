// $Id$

#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"

// constructor
NuTo::SecondPiolaKirchhoffStress1D::SecondPiolaKirchhoffStress1D()
{
    this->mSecondPiolaKirchhoffStress = 0.0;
}

// number of components
unsigned int NuTo::SecondPiolaKirchhoffStress1D::GetNumberOfComponents() const
{
    return 1;
}

// get second Piola-Kirchhoff stress
const double* NuTo::SecondPiolaKirchhoffStress1D::GetData() const
{
    return &mSecondPiolaKirchhoffStress;
}

// info routine
void NuTo::SecondPiolaKirchhoffStress1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of second Piola-Kirchhoff stress tensor (vector notation): " << this->mSecondPiolaKirchhoffStress << std::endl;
}

