// $Id$

#include <iostream>

#include "nuto/mechanics/constitutive/mechanics/RelativeInterfaceDisplacements1D.h"

NuTo::RelativeInterfaceDisplacements1D::RelativeInterfaceDisplacements1D()
{
    this->mRelativeInterfaceDisplacements = 0.0;
}

NuTo::RelativeInterfaceDisplacements1D::RelativeInterfaceDisplacements1D(const RelativeInterfaceDisplacements1D& rOther)
{
	mRelativeInterfaceDisplacements = rOther.mRelativeInterfaceDisplacements;
}

unsigned int NuTo::RelativeInterfaceDisplacements1D::GetNumberOfComponents() const
{
    return 1;
}

const double* NuTo::RelativeInterfaceDisplacements1D::GetRelativeInterfaceDisplacements1D() const
{
    return &mRelativeInterfaceDisplacements;
}

void NuTo::RelativeInterfaceDisplacements1D::GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements1D& rRelativeInterfaceDisplacements) const
{
    rRelativeInterfaceDisplacements = NuTo::RelativeInterfaceDisplacements1D(*this);
}

void NuTo::RelativeInterfaceDisplacements1D::SetRelativeInterfaceDisplacements1D(const double* rRelativeInterfaceDisplacements)
{
    this->mRelativeInterfaceDisplacements = rRelativeInterfaceDisplacements[0];
}

void NuTo::RelativeInterfaceDisplacements1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    relative interface displacement components: " << this->mRelativeInterfaceDisplacements << std::endl;
}
