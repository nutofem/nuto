// $Id$

#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/InterfaceTractions2D.h"

NuTo::InterfaceTractions2D::InterfaceTractions2D()
{
    this->mInterfaceTractions[0] = 0.0;
    this->mInterfaceTractions[1] = 0.0;
}

unsigned int NuTo::InterfaceTractions2D::GetNumberOfComponents() const
{
    return 2;
}

const double* NuTo::InterfaceTractions2D::GetData() const
{
    return this->mInterfaceTractions;
}

void NuTo::InterfaceTractions2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the interface traction vector: "
              << this->mInterfaceTractions[0] << ", "
              << this->mInterfaceTractions[1] << std::endl;
}
