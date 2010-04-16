#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"

NuTo::EngineeringStrain3D::EngineeringStrain3D()
{
	mEngineeringStrain[0] = 0.0;
	mEngineeringStrain[1] = 0.0;
	mEngineeringStrain[2] = 0.0;
	mEngineeringStrain[3] = 0.0;
	mEngineeringStrain[4] = 0.0;
	mEngineeringStrain[5] = 0.0;
}

NuTo::EngineeringStrain3D::EngineeringStrain3D(const DeformationGradient3D& rDeformationGradient)
{
    mEngineeringStrain[0] = rDeformationGradient.GetDeformationGradient3D()[0] -1;
    mEngineeringStrain[1] = rDeformationGradient.GetDeformationGradient3D()[4] -1;
    mEngineeringStrain[2] = rDeformationGradient.GetDeformationGradient3D()[8] -1;
    mEngineeringStrain[3] = rDeformationGradient.GetDeformationGradient3D()[1]+rDeformationGradient.GetDeformationGradient3D()[3];
    mEngineeringStrain[4] = rDeformationGradient.GetDeformationGradient3D()[5]+rDeformationGradient.GetDeformationGradient3D()[7];
    mEngineeringStrain[5] = rDeformationGradient.GetDeformationGradient3D()[2]+rDeformationGradient.GetDeformationGradient3D()[6];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain3D::GetNumberOfComponents() const
{
    return 6;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx,eyy,ezz,gxy,gyz,gzx)
//! @sa mDeformationGradient
const double* NuTo::EngineeringStrain3D::GetData() const
{
    return mEngineeringStrain;
}
