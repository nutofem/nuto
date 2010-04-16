#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"

NuTo::GreenLagrangeStrain3D::GreenLagrangeStrain3D()
{
	mGreenLagrangeStrain[0] = 0.0;
	mGreenLagrangeStrain[1] = 0.0;
	mGreenLagrangeStrain[2] = 0.0;
	mGreenLagrangeStrain[3] = 0.0;
	mGreenLagrangeStrain[4] = 0.0;
	mGreenLagrangeStrain[5] = 0.0;
}

NuTo::GreenLagrangeStrain3D::GreenLagrangeStrain3D(const DeformationGradient3D& rDeformationGradient)
{
	mGreenLagrangeStrain[0] = 0.5 * (rDeformationGradient.GetDeformationGradient3D()[0] * rDeformationGradient.GetDeformationGradient3D()[0] + rDeformationGradient.GetDeformationGradient3D()[1] * rDeformationGradient.GetDeformationGradient3D()[1] + rDeformationGradient.GetDeformationGradient3D()[2] * rDeformationGradient.GetDeformationGradient3D()[2] - 1.0);
	mGreenLagrangeStrain[1] = 0.5 * (rDeformationGradient.GetDeformationGradient3D()[3] * rDeformationGradient.GetDeformationGradient3D()[3] + rDeformationGradient.GetDeformationGradient3D()[4] * rDeformationGradient.GetDeformationGradient3D()[4] + rDeformationGradient.GetDeformationGradient3D()[5] * rDeformationGradient.GetDeformationGradient3D()[5] - 1.0);
	mGreenLagrangeStrain[2] = 0.5 * (rDeformationGradient.GetDeformationGradient3D()[6] * rDeformationGradient.GetDeformationGradient3D()[6] + rDeformationGradient.GetDeformationGradient3D()[7] * rDeformationGradient.GetDeformationGradient3D()[7] + rDeformationGradient.GetDeformationGradient3D()[8] * rDeformationGradient.GetDeformationGradient3D()[8] - 1.0);
	mGreenLagrangeStrain[3] = rDeformationGradient.GetDeformationGradient3D()[0] * rDeformationGradient.GetDeformationGradient3D()[3] + rDeformationGradient.GetDeformationGradient3D()[1] * rDeformationGradient.GetDeformationGradient3D()[4] + rDeformationGradient.GetDeformationGradient3D()[2] * rDeformationGradient.GetDeformationGradient3D()[5];
	mGreenLagrangeStrain[4] = rDeformationGradient.GetDeformationGradient3D()[3] * rDeformationGradient.GetDeformationGradient3D()[6] + rDeformationGradient.GetDeformationGradient3D()[4] * rDeformationGradient.GetDeformationGradient3D()[7] + rDeformationGradient.GetDeformationGradient3D()[5] * rDeformationGradient.GetDeformationGradient3D()[8];
	mGreenLagrangeStrain[5] = rDeformationGradient.GetDeformationGradient3D()[0] * rDeformationGradient.GetDeformationGradient3D()[6] + rDeformationGradient.GetDeformationGradient3D()[1] * rDeformationGradient.GetDeformationGradient3D()[7] + rDeformationGradient.GetDeformationGradient3D()[2] * rDeformationGradient.GetDeformationGradient3D()[8];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::GreenLagrangeStrain3D::GetNumberOfComponents() const
{
    return 6;
}

//! @brief ... get Green Lagrange Strain
//! @return ... Green Strain (exx,eyy,ezz,gxy,gyz,gzx)
//! @sa mDeformationGradient
const double* NuTo::GreenLagrangeStrain3D::GetData() const
{
    return mGreenLagrangeStrain;
}
