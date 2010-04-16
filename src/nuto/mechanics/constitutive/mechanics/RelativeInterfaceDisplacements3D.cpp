// $Id$

#include "nuto/mechanics/constitutive/mechanics/RelativeInterfaceDisplacements3D.h"

NuTo::RelativeInterfaceDisplacements3D::RelativeInterfaceDisplacements3D()
{
    this->mRelativeInterfaceDisplacements[0] = 0.0;
    this->mRelativeInterfaceDisplacements[1] = 0.0;
    this->mRelativeInterfaceDisplacements[2] = 0.0;
}

NuTo::RelativeInterfaceDisplacements3D::RelativeInterfaceDisplacements3D(const RelativeInterfaceDisplacements3D& rOther)
{
	mRelativeInterfaceDisplacements[0] = rOther.mRelativeInterfaceDisplacements[0];
	mRelativeInterfaceDisplacements[1] = rOther.mRelativeInterfaceDisplacements[1];
	mRelativeInterfaceDisplacements[2] = rOther.mRelativeInterfaceDisplacements[2];
}

unsigned int NuTo::RelativeInterfaceDisplacements3D::GetNumberOfComponents() const
{
    return 3;
}

const double* NuTo::RelativeInterfaceDisplacements3D::GetRelativeInterfaceDisplacements3D() const
{
    return this->mRelativeInterfaceDisplacements;
}

void NuTo::RelativeInterfaceDisplacements3D::GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements3D& rRelativeInterfaceDisplacements) const
{
	rRelativeInterfaceDisplacements = NuTo::RelativeInterfaceDisplacements3D(*this);
}

void NuTo::RelativeInterfaceDisplacements3D::SetRelativeInterfaceDisplacements3D(const double* rRelativeInterfaceDisplacements)
{
    this->mRelativeInterfaceDisplacements[0] = rRelativeInterfaceDisplacements[0];
    this->mRelativeInterfaceDisplacements[1] = rRelativeInterfaceDisplacements[1];
    this->mRelativeInterfaceDisplacements[2] = rRelativeInterfaceDisplacements[2];
}

void NuTo::RelativeInterfaceDisplacements3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    relative interface displacement components: " << this->mRelativeInterfaceDisplacements[0] << ", "
              << this->mRelativeInterfaceDisplacements[1] << ", "
              << this->mRelativeInterfaceDisplacements[2] << std::endl;
}
