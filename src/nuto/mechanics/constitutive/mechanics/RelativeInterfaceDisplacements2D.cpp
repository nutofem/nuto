// $Id$

#include "nuto/mechanics/constitutive/mechanics/RelativeInterfaceDisplacements2D.h"

NuTo::RelativeInterfaceDisplacements2D::RelativeInterfaceDisplacements2D()
{
    this->mRelativeInterfaceDisplacements[0] = 0.0;
    this->mRelativeInterfaceDisplacements[1] = 0.0;
}

NuTo::RelativeInterfaceDisplacements2D::RelativeInterfaceDisplacements2D(const RelativeInterfaceDisplacements2D& rOther)
{
	mRelativeInterfaceDisplacements[0] = rOther.mRelativeInterfaceDisplacements[0];
	mRelativeInterfaceDisplacements[1] = rOther.mRelativeInterfaceDisplacements[1];
}

unsigned int NuTo::RelativeInterfaceDisplacements2D::GetNumberOfComponents() const
{
    return 2;
}

const double* NuTo::RelativeInterfaceDisplacements2D::GetRelativeInterfaceDisplacements2D() const
{
    return this->mRelativeInterfaceDisplacements;
}

void NuTo::RelativeInterfaceDisplacements2D::GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements2D& rRelativeInterfaceDisplacements) const
{
	rRelativeInterfaceDisplacements = NuTo::RelativeInterfaceDisplacements2D(*this);
}

void NuTo::RelativeInterfaceDisplacements2D::SetRelativeInterfaceDisplacements2D(const double* rRelativeInterfaceDisplacements)
{
    this->mRelativeInterfaceDisplacements[0] = rRelativeInterfaceDisplacements[0];
    this->mRelativeInterfaceDisplacements[1] = rRelativeInterfaceDisplacements[1];
}

void NuTo::RelativeInterfaceDisplacements2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    relative interface displacement components: " << this->mRelativeInterfaceDisplacements[0] << ", "
              << this->mRelativeInterfaceDisplacements[1] << std::endl;
}
