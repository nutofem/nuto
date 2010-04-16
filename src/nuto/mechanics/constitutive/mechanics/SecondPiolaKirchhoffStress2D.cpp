// $Id$

#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"

// constructor
NuTo::SecondPiolaKirchhoffStress2D::SecondPiolaKirchhoffStress2D()
{
    this->mSecondPiolaKirchhoffStress[0] = 0.0;
    this->mSecondPiolaKirchhoffStress[1] = 0.0;
    this->mSecondPiolaKirchhoffStress[2] = 0.0;
}

// number of components
unsigned int NuTo::SecondPiolaKirchhoffStress2D::GetNumberOfComponents() const
{
    return 3;
}

// get second Piola-Kirchhoff stress
const double* NuTo::SecondPiolaKirchhoffStress2D::GetData() const
{
    return this->mSecondPiolaKirchhoffStress;
}

// info routine
void NuTo::SecondPiolaKirchhoffStress2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of second Piola-Kirchhoff stress tensor (vector notation): "
              << this->mSecondPiolaKirchhoffStress[0] << ", "
              << this->mSecondPiolaKirchhoffStress[1] << ", "
              << this->mSecondPiolaKirchhoffStress[2] << std::endl;
}

