// $Id: ConstitutiveEngineeringStressStrain.cpp 112 2009-11-17 16:43:15Z unger3 $

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"

NuTo::ConstitutiveEngineeringStressStrain::ConstitutiveEngineeringStressStrain() : ConstitutiveBase()
{
	mEnergyFlag = true;
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::ConstitutiveEngineeringStressStrain::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase);
    }
#endif // ENABLE_SERIALIZATION

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
    	throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate stress
    EngineeringStress3D stress;
    GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stress);

    // calculate strain
    EngineeringStrain3D strain;
    GetEngineeringStrain(rElement,rIp,rDeformationGradient, strain);

    //get previous data
	const ElementWithDataBase* elementWithDataBasePtr(dynamic_cast<const NuTo::ElementWithDataBase* >(rElement));
	if (elementWithDataBasePtr==0)
		throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] Element has no data.");

	const ConstitutiveStaticDataPrevEngineeringStressStrain3D*
	    staticDataPtr(dynamic_cast<const NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D* >(elementWithDataBasePtr->GetStaticData(rIp)));
	if (staticDataPtr==0)
	{
		throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] Static data is not derived from StaticDataPrevStressStrain.");
	}

	//calculate average stress
	stress.mEngineeringStress[0]=0.5*(stress.mEngineeringStress[0]+staticDataPtr->GetPrevStress()[0]);
	stress.mEngineeringStress[1]=0.5*(stress.mEngineeringStress[1]+staticDataPtr->GetPrevStress()[1]);
	stress.mEngineeringStress[2]=0.5*(stress.mEngineeringStress[2]+staticDataPtr->GetPrevStress()[2]);
	stress.mEngineeringStress[3]=0.5*(stress.mEngineeringStress[3]+staticDataPtr->GetPrevStress()[3]);
	stress.mEngineeringStress[4]=0.5*(stress.mEngineeringStress[4]+staticDataPtr->GetPrevStress()[4]);
	stress.mEngineeringStress[5]=0.5*(stress.mEngineeringStress[5]+staticDataPtr->GetPrevStress()[5]);

	//calculate delta total strain
	strain.mEngineeringStrain[0]-=staticDataPtr->GetPrevStrain()[0];
	strain.mEngineeringStrain[1]-=staticDataPtr->GetPrevStrain()[1];
	strain.mEngineeringStrain[2]-=staticDataPtr->GetPrevStrain()[2];
	strain.mEngineeringStrain[3]-=staticDataPtr->GetPrevStrain()[3];
	strain.mEngineeringStrain[4]-=staticDataPtr->GetPrevStrain()[4];
	strain.mEngineeringStrain[5]-=staticDataPtr->GetPrevStrain()[5];

	//calculate delta of the total energy
	return  staticDataPtr->GetPrevTotalEnergy()+ 0.5*(
			 stress.mEngineeringStress[0]*strain.mEngineeringStrain[0]+
	         stress.mEngineeringStress[1]*strain.mEngineeringStrain[1]+
	         stress.mEngineeringStress[2]*strain.mEngineeringStrain[2]+
	         stress.mEngineeringStress[3]*strain.mEngineeringStrain[3]+
	         stress.mEngineeringStress[4]*strain.mEngineeringStrain[4]+
	         stress.mEngineeringStress[5]*strain.mEngineeringStrain[5]);
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
	// check if parameters are valid
	if (this->mParametersValid == false)
	{
		throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
	}

	// calculate stress
	EngineeringStress3D stress;
	GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stress);

	// calculate delta strain
	EngineeringStrain3D deltaElasticStrain;
	GetDeltaElasticEngineeringStrain(rElement,rIp,rDeformationGradient, deltaElasticStrain);

	//get previous data
	const ElementWithDataBase* elementWithDataBasePtr(dynamic_cast<const NuTo::ElementWithDataBase* >(rElement));
	if (elementWithDataBasePtr==0)
		throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] Element has no data.");

	const ConstitutiveStaticDataPrevEngineeringStressStrain3D*
		staticDataPtr(dynamic_cast<const NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D* >(elementWithDataBasePtr->GetStaticData(rIp)));
	if (staticDataPtr==0)
	{
		throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain] Static data is not derived from StaticDataPrevStressStrain.");
	}

	//calculate average stress
	stress.mEngineeringStress[0]=0.5*(stress.mEngineeringStress[0]+staticDataPtr->GetPrevStress()[0]);
	stress.mEngineeringStress[1]=0.5*(stress.mEngineeringStress[1]+staticDataPtr->GetPrevStress()[1]);
	stress.mEngineeringStress[2]=0.5*(stress.mEngineeringStress[2]+staticDataPtr->GetPrevStress()[2]);
	stress.mEngineeringStress[3]=0.5*(stress.mEngineeringStress[3]+staticDataPtr->GetPrevStress()[3]);
	stress.mEngineeringStress[4]=0.5*(stress.mEngineeringStress[4]+staticDataPtr->GetPrevStress()[4]);
	stress.mEngineeringStress[5]=0.5*(stress.mEngineeringStress[5]+staticDataPtr->GetPrevStress()[5]);

	//calculate delta of the elastic energy
	return  staticDataPtr->GetPrevElasticEnergy()+ 0.5*(
			 stress.mEngineeringStress[0]*deltaElasticStrain.mEngineeringStrain[0]+
			 stress.mEngineeringStress[1]*deltaElasticStrain.mEngineeringStrain[1]+
			 stress.mEngineeringStress[2]*deltaElasticStrain.mEngineeringStrain[2]+
			 stress.mEngineeringStress[3]*deltaElasticStrain.mEngineeringStrain[3]+
			 stress.mEngineeringStress[4]*deltaElasticStrain.mEngineeringStrain[4]+
			 stress.mEngineeringStress[5]*deltaElasticStrain.mEngineeringStrain[5]);
}
