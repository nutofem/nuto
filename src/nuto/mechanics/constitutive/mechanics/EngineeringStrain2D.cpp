#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"

NuTo::EngineeringStrain2D::EngineeringStrain2D()
{
	mEngineeringStrain[0] = 0.0;
	mEngineeringStrain[1] = 0.0;
	mEngineeringStrain[2] = 0.0;
}

NuTo::EngineeringStrain2D::EngineeringStrain2D(const DeformationGradient2D& rDeformationGradient)
{
    mEngineeringStrain[0] = rDeformationGradient.GetDeformationGradient2D()[0] -1;
    mEngineeringStrain[1] = rDeformationGradient.GetDeformationGradient2D()[3] -1;;
    mEngineeringStrain[2] = rDeformationGradient.GetDeformationGradient2D()[1]+rDeformationGradient.GetDeformationGradient2D()[2];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain2D::GetNumberOfComponents() const
{
    return 3;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx,eyy,gxy)
//! @sa mDeformationGradient
const double* NuTo::EngineeringStrain2D::GetData() const
{
    return mEngineeringStrain;
}
