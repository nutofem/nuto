#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"

NuTo::GreenLagrangeStrain1D::GreenLagrangeStrain1D()
{
	mGreenLagrangeStrain = 0.0;
}

NuTo::GreenLagrangeStrain1D::GreenLagrangeStrain1D(const DeformationGradient1D& rDeformationGradient)
{
    mGreenLagrangeStrain = 0.5*(rDeformationGradient.GetDeformationGradient1D()[0]*rDeformationGradient.GetDeformationGradient1D()[0] -1);
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::GreenLagrangeStrain1D::GetNumberOfComponents() const
{
    return 1;
}

//! @brief ... get Green Strain
//! @return ... Green Strain (exx)
const double* NuTo::GreenLagrangeStrain1D::GetData() const
{
    return &mGreenLagrangeStrain;
}
