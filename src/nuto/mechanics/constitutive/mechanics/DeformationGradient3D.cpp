// $Id$

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"

// constructor
NuTo::DeformationGradient3D::DeformationGradient3D()
{
    mDeformationGradient[0] = 0.0;
    mDeformationGradient[1] = 0.0;
    mDeformationGradient[2] = 0.0;
    mDeformationGradient[3] = 0.0;
    mDeformationGradient[4] = 0.0;
    mDeformationGradient[5] = 0.0;
    mDeformationGradient[6] = 0.0;
    mDeformationGradient[7] = 0.0;
    mDeformationGradient[8] = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient3D::DeformationGradient3D(const DeformationGradient1D& rOther)
{
	mDeformationGradient[0] = rOther.GetDeformationGradient1D()[0];
    mDeformationGradient[1] = 0.0;
    mDeformationGradient[2] = 0.0;
    mDeformationGradient[3] = 0.0;
    mDeformationGradient[4] = 0.0;
    mDeformationGradient[5] = 0.0;
    mDeformationGradient[6] = 0.0;
    mDeformationGradient[7] = 0.0;
    mDeformationGradient[8] = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient3D::DeformationGradient3D(const DeformationGradient2D& rOther)
{
    mDeformationGradient[0] = rOther.GetDeformationGradient2D()[0];
    mDeformationGradient[1] = rOther.GetDeformationGradient2D()[1];
    mDeformationGradient[2] = 0.0;
    mDeformationGradient[3] = rOther.GetDeformationGradient2D()[2];
    mDeformationGradient[4] = rOther.GetDeformationGradient2D()[3];
    mDeformationGradient[5] = 0.0;
    mDeformationGradient[6] = 0.0;
    mDeformationGradient[7] = 0.0;
    mDeformationGradient[8] = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient3D::DeformationGradient3D(const DeformationGradient3D& rOther)
{
	mDeformationGradient[0] = rOther.mDeformationGradient[0];
    mDeformationGradient[1] = rOther.mDeformationGradient[1];
    mDeformationGradient[2] = rOther.mDeformationGradient[2];
    mDeformationGradient[3] = rOther.mDeformationGradient[3];
    mDeformationGradient[4] = rOther.mDeformationGradient[4];
    mDeformationGradient[5] = rOther.mDeformationGradient[5];
    mDeformationGradient[6] = rOther.mDeformationGradient[6];
    mDeformationGradient[7] = rOther.mDeformationGradient[7];
    mDeformationGradient[8] = rOther.mDeformationGradient[8];
}

// get number of components
unsigned int NuTo::DeformationGradient3D::GetNumberOfComponents() const
{
    return 9;
}

// get deformation gradient
const double* NuTo::DeformationGradient3D::GetDeformationGradient3D() const
{
    return this->mDeformationGradient;
}

void NuTo::DeformationGradient3D::GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const
{
	rDeformationGradient = NuTo::DeformationGradient3D(*this);
}
// set deformation gradient
void NuTo::DeformationGradient3D::SetDeformationGradient3D(const double* rDeformationGradient)
{
    for (unsigned int count = 0; count < 9; count++)
    {
        this->mDeformationGradient[count] = rDeformationGradient[count];
    }
}

// calculate engineering strain
void NuTo::DeformationGradient3D::GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const
{
	rEngineeringStrain = NuTo::EngineeringStrain3D(*this);
}

// calculate Green strain
void NuTo::DeformationGradient3D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenLagrangeStrain) const
{
	rGreenLagrangeStrain = NuTo::GreenLagrangeStrain3D(*this);
}
