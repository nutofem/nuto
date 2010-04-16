// $Id$

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"

#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"

#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"

// constructor
NuTo::DeformationGradient1D::DeformationGradient1D()
{
    this->mDeformationGradient = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient1D::DeformationGradient1D(const DeformationGradient1D& rOther)
{
	mDeformationGradient = rOther.mDeformationGradient;
}

// get number of components
unsigned int NuTo::DeformationGradient1D::GetNumberOfComponents() const
{
    return 1;
}

// get deformation gradient
const double* NuTo::DeformationGradient1D::GetDeformationGradient1D() const
{
    return &mDeformationGradient;
}

// get deformation gradient
void NuTo::DeformationGradient1D::GetDeformationGradient(NuTo::DeformationGradient1D& rDeformationGradient) const
{
    rDeformationGradient = NuTo::DeformationGradient1D(*this);
}

// get deformation gradient
void NuTo::DeformationGradient1D::GetDeformationGradient(NuTo::DeformationGradient2D& rDeformationGradient) const
{
    rDeformationGradient = NuTo::DeformationGradient2D(*this);
}


void NuTo::DeformationGradient1D::GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const
{
	rDeformationGradient = NuTo::DeformationGradient3D(*this);}

// set deformation gradient
void NuTo::DeformationGradient1D::SetDeformationGradient1D(const double* rDeformationGradient)
{
    this->mDeformationGradient = rDeformationGradient[0];
}

// calculate engineering strain
void NuTo::DeformationGradient1D::GetEngineeringStrain(NuTo::EngineeringStrain1D& rEngineeringStrain) const
{
    rEngineeringStrain = NuTo::EngineeringStrain1D(*this);
}

void NuTo::DeformationGradient1D::GetEngineeringStrain(NuTo::EngineeringStrain2D& rEngineeringStrain) const
{
	rEngineeringStrain = NuTo::EngineeringStrain2D(*this);
}

void NuTo::DeformationGradient1D::GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const
{
	rEngineeringStrain = NuTo::EngineeringStrain3D(*this);
}

// calculate Green strain
void NuTo::DeformationGradient1D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain1D& rGreenLagrangeStrain) const
{
    rGreenLagrangeStrain = NuTo::GreenLagrangeStrain1D(*this);
}

void NuTo::DeformationGradient1D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain2D& rGreenLagrangeStrain) const
{
	rGreenLagrangeStrain = NuTo::GreenLagrangeStrain2D(*this);
}

void NuTo::DeformationGradient1D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenLagrangeStrain) const
{
	rGreenLagrangeStrain = NuTo::GreenLagrangeStrain3D(*this);
}
