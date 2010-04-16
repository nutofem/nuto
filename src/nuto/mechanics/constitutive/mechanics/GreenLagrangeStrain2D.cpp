#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"

NuTo::GreenLagrangeStrain2D::GreenLagrangeStrain2D()
{
	mGreenLagrangeStrain[0] = 0.0;
	mGreenLagrangeStrain[1] = 0.0;
	mGreenLagrangeStrain[2] = 0.0;
}


NuTo::GreenLagrangeStrain2D::GreenLagrangeStrain2D(const DeformationGradient2D& rDeformationGradient)
{
    mGreenLagrangeStrain[0] = 0.5 * (rDeformationGradient.GetDeformationGradient2D()[0] * rDeformationGradient.GetDeformationGradient2D()[0] + rDeformationGradient.GetDeformationGradient2D()[1] * rDeformationGradient.GetDeformationGradient2D()[1] - 1.0);
    mGreenLagrangeStrain[1] = 0.5 * (rDeformationGradient.GetDeformationGradient2D()[3] * rDeformationGradient.GetDeformationGradient2D()[3] + rDeformationGradient.GetDeformationGradient2D()[2] * rDeformationGradient.GetDeformationGradient2D()[2] - 1.0);
    mGreenLagrangeStrain[2] = rDeformationGradient.GetDeformationGradient2D()[0] * rDeformationGradient.GetDeformationGradient2D()[2] + rDeformationGradient.GetDeformationGradient2D()[1] * rDeformationGradient.GetDeformationGradient2D()[3];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::GreenLagrangeStrain2D::GetNumberOfComponents() const
{
    return 3;
}

//! @brief ... get green lagrange Strain
//! @return ... (exx,eyy,gxy)
//! @sa mDeformationGradient
const double* NuTo::GreenLagrangeStrain2D::GetData() const
{
    return mGreenLagrangeStrain;
}
