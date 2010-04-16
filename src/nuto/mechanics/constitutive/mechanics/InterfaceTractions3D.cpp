// $Id$

#include <iostream>
#include "nuto/mechanics/constitutive/mechanics/InterfaceTractions3D.h"

NuTo::InterfaceTractions3D::InterfaceTractions3D()
{
    this->mInterfaceTractions[0] = 0.0;
    this->mInterfaceTractions[1] = 0.0;
    this->mInterfaceTractions[2] = 0.0;
}

unsigned int NuTo::InterfaceTractions3D::GetNumberOfComponents() const
{
    return 3;
}

const double* NuTo::InterfaceTractions3D::GetData() const
{
    return this->mInterfaceTractions;
}

void NuTo::InterfaceTractions3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of the interface traction vector: "
              << this->mInterfaceTractions[0] << ", "
              << this->mInterfaceTractions[1] << ", "
              << this->mInterfaceTractions[2] << std::endl;
}
