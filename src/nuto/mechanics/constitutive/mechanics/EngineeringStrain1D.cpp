#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"

NuTo::EngineeringStrain1D::EngineeringStrain1D()
{
	mEngineeringStrain = 0.0;
}

NuTo::EngineeringStrain1D::EngineeringStrain1D(const DeformationGradient1D& rDeformationGradient)
{
    mEngineeringStrain = rDeformationGradient.GetDeformationGradient1D()[0] -1;
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain1D::GetNumberOfComponents() const
{
    return 1;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx)
const double* NuTo::EngineeringStrain1D::GetData() const
{
    return &mEngineeringStrain;
}
