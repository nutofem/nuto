// $Id: ConstitutiveEngineeringStressStrain.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::ConstitutiveEngineeringStressStrain::ConstitutiveEngineeringStressStrain() : ConstitutiveBase()
{
	mEnergyFlag = true;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveEngineeringStressStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveEngineeringStressStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveEngineeringStressStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveEngineeringStressStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveEngineeringStressStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveEngineeringStressStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveEngineeringStressStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveEngineeringStressStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
       & BOOST_SERIALIZATION_NVP(mEnergyFlag);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveEngineeringStressStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveEngineeringStressStrain)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstitutiveEngineeringStressStrain)
#endif // ENABLE_SERIALIZATION

    //  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering strain from deformation gradient in 1D (truss is assumed to be plane stress)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const
{
	EngineeringStrain1D engineeringStrain;
	rDeformationGradient.GetEngineeringStrain(engineeringStrain);

	double mNu(GetPoissonsRatio());

	rEngineeringStrain.mEngineeringStrain[0] = engineeringStrain.mEngineeringStrain;
	rEngineeringStrain.mEngineeringStrain[1] = -mNu*engineeringStrain.mEngineeringStrain;
	rEngineeringStrain.mEngineeringStrain[2] = -mNu*engineeringStrain.mEngineeringStrain;
	rEngineeringStrain.mEngineeringStrain[3] = 0.;
	rEngineeringStrain.mEngineeringStrain[4] = 0.;
	rEngineeringStrain.mEngineeringStrain[5] = 0.;

	return Error::SUCCESSFUL;
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering strain from deformation gradient in 3D
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetEngineeringStrain(const ElementBase* rElement,const int rIp,
		      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const
{
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    assert(rElement->GetSection()!=0);

    switch(rElement->GetSection()->GetType())
    {
    case Section::PLANE_STRAIN:
		rEngineeringStrain.mEngineeringStrain[0] = engineeringStrain.mEngineeringStrain[0];
		rEngineeringStrain.mEngineeringStrain[1] = engineeringStrain.mEngineeringStrain[1];
		rEngineeringStrain.mEngineeringStrain[2] = 0;
		rEngineeringStrain.mEngineeringStrain[3] = engineeringStrain.mEngineeringStrain[2];
		rEngineeringStrain.mEngineeringStrain[4] = 0.;
		rEngineeringStrain.mEngineeringStrain[5] = 0.;
    	break;
    case Section::PLANE_STRESS:
		GetEngineeringStrainFromEngineeringStrain(rElement, rIp,engineeringStrain,rEngineeringStrain);
    	break;
    default:
    	throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetEngineeringStrain] Wrong type of 2D section type");
    }

	return Error::SUCCESSFUL;
}


//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering strain from deformation gradient in 3D
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const
{
	rDeformationGradient.GetEngineeringStrain(rEngineeringStrain);
	return Error::SUCCESSFUL;
}

//  Engineering strain2D -> 3D /////////////////////////////////////
//! @brief ... calculate engineering strain from deformation gradient in 3D
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetEngineeringStrainFromEngineeringStrain(const ElementBase* rElement, int rIp,
								  const EngineeringStrain2D& rEngineeringStrain2D, EngineeringStrain3D& rEngineeringStrain3D) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetEngineeringStrainFromEngineeringStrain] Not implemented for this material law");
}

//  Energy /////////////////////////////////////
//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetInternalEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetInternalEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient, double& rEnergy) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
    	throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::GetInternalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate stress
    EngineeringStress3D stress;
    GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stress);

    // calculate strain
    EngineeringStrain3D strain;
    GetEngineeringStrain(rElement,rIp,rDeformationGradient, strain);

    //get previous data
	const ConstitutiveStaticDataPrevEngineeringStressStrain3D*
	    staticDataPtr(dynamic_cast<const NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D* >(rElement->GetStaticData(rIp)));
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
	rEnergy = staticDataPtr->GetPrevTotalEnergy()+ 0.5*(
			 stress.mEngineeringStress[0]*strain.mEngineeringStrain[0]+
	         stress.mEngineeringStress[1]*strain.mEngineeringStrain[1]+
	         stress.mEngineeringStress[2]*strain.mEngineeringStrain[2]+
	         stress.mEngineeringStress[3]*strain.mEngineeringStrain[3]+
	         stress.mEngineeringStress[4]*strain.mEngineeringStrain[4]+
	         stress.mEngineeringStress[5]*strain.mEngineeringStrain[5]);
    return Error::SUCCESSFUL;
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain] to be implemented.");
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, double& rEnergy) const
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
	const ConstitutiveStaticDataPrevEngineeringStressStrain3D*
		staticDataPtr(dynamic_cast<const NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D* >(rElement->GetStaticData(rIp)));
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
	rEnergy = staticDataPtr->GetPrevElasticEnergy()+ 0.5*(
			 stress.mEngineeringStress[0]*deltaElasticStrain.mEngineeringStrain[0]+
			 stress.mEngineeringStress[1]*deltaElasticStrain.mEngineeringStrain[1]+
			 stress.mEngineeringStress[2]*deltaElasticStrain.mEngineeringStrain[2]+
			 stress.mEngineeringStress[3]*deltaElasticStrain.mEngineeringStrain[3]+
			 stress.mEngineeringStress[4]*deltaElasticStrain.mEngineeringStrain[4]+
			 stress.mEngineeringStress[5]*deltaElasticStrain.mEngineeringStrain[5]);
    return Error::SUCCESSFUL;
}
//! @brief ... avoid dynamic cast
//! @return ... see brief explanation
NuTo::ConstitutiveEngineeringStressStrain* NuTo::ConstitutiveEngineeringStressStrain::AsConstitutiveEngineeringStressStrain()
{
	return this;
}

//! @brief ... avoid dynamic cast
//! @return ... see brief explanation
const NuTo::ConstitutiveEngineeringStressStrain* NuTo::ConstitutiveEngineeringStressStrain::AsConstitutiveEngineeringStressStrain()const
{
	return this;
}

//! @brief ... checks, if a model has to be switched from linear to nonlinear, and then performs the adaption
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveEngineeringStressStrain::MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient)const
{
	throw MechanicsException("[NuTo::ConstitutiveEngineeringStressStrain::MultiscaleSwitchToNonlinear] not implemented for this material routine.");
}
