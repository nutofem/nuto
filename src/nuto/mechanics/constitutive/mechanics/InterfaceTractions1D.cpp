// $Id$

#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/InterfaceTractions1D.h"

NuTo::InterfaceTractions1D::InterfaceTractions1D()
{
    this->mInterfaceTractions = 0.0;
}

unsigned int NuTo::InterfaceTractions1D::GetNumberOfComponents() const
{
    return 1;
}

const double* NuTo::InterfaceTractions1D::GetData() const
{
    return &mInterfaceTractions;
}


void NuTo::InterfaceTractions1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the interface traction vector: " << this->mInterfaceTractions << std::endl;
}
