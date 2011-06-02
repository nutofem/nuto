// $Id: ConstitutiveMisesPlasticity.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveMisesPlasticity.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticity3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticityWithEnergy3D.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/math/FullMatrix.h"

NuTo::ConstitutiveMisesPlasticity::ConstitutiveMisesPlasticity() : ConstitutiveEngineeringStressStrain()
{
	mE = 0.;
	mNu = 0.;
	mSigma.resize(1);
	mH.resize(1);
	SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveMisesPlasticity::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveMisesPlasticity::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveMisesPlasticity::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveMisesPlasticity::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveMisesPlasticity::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveMisesPlasticity::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveMisesPlasticity::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    rLogger << "start serialize ConstitutiveMisesPlasticity" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
       & BOOST_SERIALIZATION_NVP(mE)
       & BOOST_SERIALIZATION_NVP(mNu)
       & BOOST_SERIALIZATION_NVP(mSigma)
       & BOOST_SERIALIZATION_NVP(mH);
#ifdef DEBUG_SERIALIZATION
    rLogger << "finish serialize ConstitutiveMisesPlasticity" << "\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveMisesPlasticity)
#endif // ENABLE_SERIALIZATION

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStrain] To be implemented.");
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStrain] To be implemented.");
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    // perform return mapping
	ConstitutiveStaticDataMisesPlasticity3D  staticData;
    ReturnMapping3D(rElement, rIp, rDeformationGradient, 0 , 0, &staticData,rElement->GetStructure()->GetLogger());
    rEngineeringPlasticStrain.mEngineeringStrain[0] = staticData.mEpsilonP[0];
    rEngineeringPlasticStrain.mEngineeringStrain[1] = staticData.mEpsilonP[1];
    rEngineeringPlasticStrain.mEngineeringStrain[2] = staticData.mEpsilonP[2];
    rEngineeringPlasticStrain.mEngineeringStrain[3] = staticData.mEpsilonP[3];
    rEngineeringPlasticStrain.mEngineeringStrain[4] = staticData.mEpsilonP[4];
    rEngineeringPlasticStrain.mEngineeringStrain[5] = staticData.mEpsilonP[5];
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStrain] To be implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStrain] To be implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain] To be implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain] To be implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain3D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // perform return mapping
    ReturnMapping3D(rElement, rIp, rDeformationGradient, &rEngineeringStress, 0, 0,rElement->GetStructure()->GetLogger());
}

//  Damage /////////////////////////////////////
 //! @brief ... calculate isotropic damage from deformation gradient in 1D
 //! @param rElement ... element
 //! @param rIp ... integration point
 //! @param rDeformationGradient ... deformation gradient
 //! @param rDamage ... damage variable
 void NuTo::ConstitutiveMisesPlasticity::GetDamage(const ElementBase* rElement, int rIp,
                                   const DeformationGradient1D& rDeformationGradient, double& rDamage) const
{
	rDamage=0.;
}

 //  Damage /////////////////////////////////////
 //! @brief ... calculate isotropic damage from deformation gradient in 2D
 //! @param rElement ... element
 //! @param rIp ... integration point
 //! @param rDeformationGradient ... deformation gradient
 //! @param rDamage ... damage variable
 void NuTo::ConstitutiveMisesPlasticity::GetDamage(const ElementBase* rElement, int rIp,
                                   const DeformationGradient2D& rDeformationGradient, double& rDamage) const
{
    rDamage=0.;
}

 //  Damage /////////////////////////////////////
 //! @brief ... calculate isotropic damage from deformation gradient in 3D
 //! @param rElement ... element
 //! @param rIp ... integration point
 //! @param rDeformationGradient ... deformation gradient
 //! @param rDamage ... damage variable
 void NuTo::ConstitutiveMisesPlasticity::GetDamage(const ElementBase* rElement, int rIp,
                                   const DeformationGradient3D& rDeformationGradient, double& rDamage) const
{
    rDamage=0.;
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::ConstitutiveMisesPlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient,
		ConstitutiveTangentBase* rTangent) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetTangent_EngineeringStress_EngineeringStrain] To be implemented.");
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::ConstitutiveMisesPlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient,
		ConstitutiveTangentBase* rTangent) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetTangent_EngineeringStress_EngineeringStrain] To be implemented.");
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::ConstitutiveMisesPlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient,
		ConstitutiveTangentBase* rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // perform return mapping
    ReturnMapping3D(rElement, rIp, rDeformationGradient, 0 , rTangent->AsConstitutiveTangentLocal6x6(), 0,rElement->GetStructure()->GetLogger());
    rTangent->SetSymmetry(true);
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::ConstitutiveMisesPlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient)const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain] To be implemented.");
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::ConstitutiveMisesPlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient)const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain] To be implemented.");
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::ConstitutiveMisesPlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient)const
{
    // perform return mapping
	ConstitutiveStaticDataMisesPlasticity3D*
	    staticDataPtr(dynamic_cast<NuTo::ConstitutiveStaticDataMisesPlasticity3D* >(rElement->GetStaticData(rIp)));
	if (staticDataPtr==0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain] Static data is not derived from Mises3D.");
    ReturnMapping3D(rElement, rIp, rDeformationGradient, 0 , 0, staticDataPtr,rElement->GetStructure()->GetLogger());
    rElement->GetStructure()->GetLogger()<<"updated plastic strain " << staticDataPtr->mEpsilonP[0] << "\n";
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::ConstitutiveMisesPlasticity::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient)const
{
	//no need to update tmp static data
}


//! @brief ... update mp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::ConstitutiveMisesPlasticity::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient)const
{
	//no need to update tmp static data
}

//! @brief ... update mp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::ConstitutiveMisesPlasticity::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient)const
{
	//no need to update tmp static data
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveMisesPlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
		const ElementBase* rElement) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D] To be implemented.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveMisesPlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
		const ElementBase* rElement) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain2D] To be implemented.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveMisesPlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
		const ElementBase* rElement) const
{
	if (mEnergyFlag)
	{
		return new NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D();
	}
	else
	{
		return new NuTo::ConstitutiveStaticDataMisesPlasticity3D();
	}
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::ConstitutiveMisesPlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetDeltaElasticEngineeringStrain] To be implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::ConstitutiveMisesPlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetDeltaElasticEngineeringStrain] To be implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::ConstitutiveMisesPlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
    // calculate current strain
    EngineeringStrain3D curTotalEngineeringStrain;
    rDeformationGradient.GetEngineeringStrain(curTotalEngineeringStrain);

	// perform return mapping to calculate current plastic strain
	ConstitutiveStaticDataMisesPlasticity3D newStaticData;
    ReturnMapping3D(rElement, rIp, rDeformationGradient, 0 , 0, &newStaticData,rElement->GetStructure()->GetLogger());

	//get previous data
	const ConstitutiveStaticDataMisesPlasticityWithEnergy3D*
		staticDataPtr(dynamic_cast<const NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D* >(rElement->GetStaticData(rIp)));
	if (staticDataPtr==0)
	{
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetDeltaElasticEngineeringStrain] Static data is not derived from ConstitutiveStaticDataMisesPlasticityWithEnergy3D.");
	}

	//current elastic strain (which finally given delta_elastic = cur_elastic-prev_elastic)
	rDeltaElasticEngineeringStrain.mEngineeringStrain[0] = curTotalEngineeringStrain.mEngineeringStrain[0] - newStaticData.mEpsilonP[0];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[1] = curTotalEngineeringStrain.mEngineeringStrain[1] - newStaticData.mEpsilonP[1];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[2] = curTotalEngineeringStrain.mEngineeringStrain[2] - newStaticData.mEpsilonP[2];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[3] = curTotalEngineeringStrain.mEngineeringStrain[3] - newStaticData.mEpsilonP[3];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[4] = curTotalEngineeringStrain.mEngineeringStrain[4] - newStaticData.mEpsilonP[4];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[5] = curTotalEngineeringStrain.mEngineeringStrain[5] - newStaticData.mEpsilonP[5];

	//calculate previous elastic strain with a minus sign
	rDeltaElasticEngineeringStrain.mEngineeringStrain[0] -= staticDataPtr->GetPrevStrain()[0] - staticDataPtr->mEpsilonP[0];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[1] -= staticDataPtr->GetPrevStrain()[1] - staticDataPtr->mEpsilonP[1];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[2] -= staticDataPtr->GetPrevStrain()[2] - staticDataPtr->mEpsilonP[2];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[3] -= staticDataPtr->GetPrevStrain()[3] - staticDataPtr->mEpsilonP[3];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[4] -= staticDataPtr->GetPrevStrain()[4] - staticDataPtr->mEpsilonP[4];
	rDeltaElasticEngineeringStrain.mEngineeringStrain[5] -= staticDataPtr->GetPrevStrain()[5] - staticDataPtr->mEpsilonP[5];
}

//! @brief ... performs the return mapping procedure in 3D
//! @param rElement ... structure
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient point
//! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
//! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
//! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
#define sqrt_2div3 0.81649658
#define tolerance 1e-8
void NuTo::ConstitutiveMisesPlasticity::ReturnMapping3D(const ElementBase* rElement,int rIp,
		const DeformationGradient3D& rDeformationGradient,
		EngineeringStress3D* rNewStress,
		ConstitutiveTangentLocal6x6* rNewTangent,
		ConstitutiveStaticDataMisesPlasticity3D* rNewStaticData,
		Logger& rLogger)const
{
    double sigma_trial[6],
    xi_trial[6],
    norm_dev,
    sigma_y,
    factor,
    factor2,
    yield_condition,
    epsilon_p_eq2,
    d_sigma,
    d_H,
    H,
    H2,
    mu,
    bulk_modulus,
    delta_gamma=0.,
    g,
    dg,
    df_dsigma[6],
    trace_epsilon,
    trace_epsilon_div_3;

    //here use the modified parameters to enable random fields
    double modE = GetRanfieldFactorYoungsModulus(rElement,rIp)*mE;
    CheckYoungsModulus(modE);

    double modNu = GetRanfieldFactorPoissonsRatio(rElement,rIp)*mNu;
    CheckPoissonsRatio(modNu);

    double ranfield_factor_yield_strength(GetRanfieldFactorYieldStrength(rElement,rIp));
    double ranfield_factor_hardening_modulus(GetRanfieldFactorHardeningModulus(rElement,rIp));

    mu    = modE/(2.*(1.+modNu));
    bulk_modulus = modE/(3.-6.*modNu);

    // calculate current strain
    EngineeringStrain3D curTotalEngineeringStrain;
    rDeformationGradient.GetEngineeringStrain(curTotalEngineeringStrain);
    const double *total_strain(&(curTotalEngineeringStrain.mEngineeringStrain[0]));

    //get old static data
    const ConstitutiveStaticDataMisesPlasticity3D* rOldStaticData;
    rOldStaticData = dynamic_cast<const ConstitutiveStaticDataMisesPlasticity3D*>(rElement->GetStaticData(rIp));
    if (rOldStaticData==0)
    	throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::ReturnMapping3D] Static data pointer is not of type ConstitutiveStaticDataMisesPlasticity3D.");
    const double *plastic_strain(rOldStaticData->mEpsilonP);
    const double *back_stress(rOldStaticData->mSigmaB);

    trace_epsilon = total_strain[0] + total_strain[1] + total_strain[2];

    trace_epsilon_div_3 = trace_epsilon/3.;

    //trial stress
    sigma_trial[0] = (total_strain[0]-trace_epsilon_div_3-plastic_strain[0])*2.*mu;
    sigma_trial[1] = (total_strain[1]-trace_epsilon_div_3-plastic_strain[1])*2.*mu;
    sigma_trial[2] = (total_strain[2]-trace_epsilon_div_3-plastic_strain[2])*2.*mu;
    sigma_trial[3] = (total_strain[3]                    -plastic_strain[3])*mu; // in total strain, gamma is stored
    sigma_trial[4] = (total_strain[4]                    -plastic_strain[4])*mu; // in total strain, gamma is stored
    sigma_trial[5] = (total_strain[5]                    -plastic_strain[5])*mu; // in total strain, gamma is stored

    //subtract backstress
    xi_trial[0] = sigma_trial[0]-back_stress[0];
    xi_trial[1] = sigma_trial[1]-back_stress[1];
    xi_trial[2] = sigma_trial[2]-back_stress[2];
    xi_trial[3] = sigma_trial[3]-back_stress[3];
    xi_trial[4] = sigma_trial[4]-back_stress[4];
    xi_trial[5] = sigma_trial[5]-back_stress[5];

    //printf("total_strain %g %g %g %g %g %g\n",total_strain[0],total_strain[1],total_strain[2],total_strain[3],total_strain[4],total_strain[5]);
    //printf("Xi trial %g %g %g %g %g %g\n",xi_trial[0],xi_trial[1],xi_trial[2],xi_trial[3],xi_trial[4],xi_trial[5]);

    // norm of deviator
    norm_dev = sqrt(xi_trial[0]*xi_trial[0]+xi_trial[1]*xi_trial[1]+xi_trial[2]*xi_trial[2]+
                    2.*(xi_trial[3]*xi_trial[3]+xi_trial[4]*xi_trial[4]+xi_trial[5]*xi_trial[5]));

    //determine radius of yield function
    sigma_y = ranfield_factor_yield_strength*GetYieldStrength(rOldStaticData->mEpsilonPEq,d_sigma);
    d_sigma*= ranfield_factor_yield_strength;
    yield_condition = norm_dev - sqrt_2div3 * sigma_y;

    if (yield_condition<-tolerance*sigma_y)
    {
        // elastic regime
    	factor = bulk_modulus*trace_epsilon;
        if (rNewStress!=0)
        {
        	rNewStress->mEngineeringStress[0] = factor+sigma_trial[0];
        	rNewStress->mEngineeringStress[1] = factor+sigma_trial[1];
        	rNewStress->mEngineeringStress[2] = factor+sigma_trial[2];
        	rNewStress->mEngineeringStress[3] = 	   sigma_trial[3];
        	rNewStress->mEngineeringStress[4] = 	   sigma_trial[4];
        	rNewStress->mEngineeringStress[5] = 	   sigma_trial[5];
        }
        if (rNewTangent!=0)
        {
            factor = modE/(1.+modNu)/(1.-2.*modNu);

            rNewTangent->mTangent[0] = (1.-modNu)*factor;
            rNewTangent->mTangent[1] = modNu*factor;
            rNewTangent->mTangent[2] = rNewTangent->mTangent[1];
            rNewTangent->mTangent[3] = 0.;
            rNewTangent->mTangent[4] = 0.;
            rNewTangent->mTangent[5] = 0.;

            rNewTangent->mTangent[6] = rNewTangent->mTangent[1];
            rNewTangent->mTangent[7] = rNewTangent->mTangent[0];
            rNewTangent->mTangent[8] = rNewTangent->mTangent[1];
            rNewTangent->mTangent[9] = 0.;
            rNewTangent->mTangent[10] = 0.;
            rNewTangent->mTangent[11] = 0.;

            rNewTangent->mTangent[12] = rNewTangent->mTangent[1];
            rNewTangent->mTangent[13] = rNewTangent->mTangent[1];
            rNewTangent->mTangent[14] = rNewTangent->mTangent[0];
            rNewTangent->mTangent[15] = 0.;
            rNewTangent->mTangent[16] = 0.;
            rNewTangent->mTangent[17] = 0.;

            rNewTangent->mTangent[18] = 0.;
            rNewTangent->mTangent[19] = 0.;
            rNewTangent->mTangent[20] = 0.;
            rNewTangent->mTangent[21] = (0.5-modNu)*factor;
            rNewTangent->mTangent[22] = 0.;
            rNewTangent->mTangent[23] = 0.;

            rNewTangent->mTangent[24] = 0.;
            rNewTangent->mTangent[25] = 0.;
            rNewTangent->mTangent[26] = 0.;
            rNewTangent->mTangent[27] = 0.;
            rNewTangent->mTangent[28] = rNewTangent->mTangent[21];
            rNewTangent->mTangent[29] = 0.;

            rNewTangent->mTangent[30] = 0.;
            rNewTangent->mTangent[31] = 0.;
            rNewTangent->mTangent[32] = 0.;
            rNewTangent->mTangent[33] = 0.;
            rNewTangent->mTangent[34] = 0.;
            rNewTangent->mTangent[35] = rNewTangent->mTangent[21];
        }
        // static data is unchanged
        return;
    }

    //plastic loading
    H  = ranfield_factor_hardening_modulus * GetHardeningModulus(rOldStaticData->mEpsilonPEq,d_H);
    d_H *= ranfield_factor_hardening_modulus;
    epsilon_p_eq2 = rOldStaticData->mEpsilonPEq;
    H2 = H;

    int i=0;
    for (;i<100;i++)
    {
        g  = yield_condition - (2.*mu*delta_gamma + sqrt_2div3 * (H2-H));
        if (fabs(g)<tolerance*sigma_y)
        {
            break;
        }
        dg = -2.* mu * (1.+(d_H+d_sigma)/(3.*mu));
        delta_gamma -=g/dg;
        epsilon_p_eq2 = rOldStaticData->mEpsilonPEq + sqrt_2div3 * delta_gamma;

        H2  = ranfield_factor_hardening_modulus * GetHardeningModulus(epsilon_p_eq2,d_H);
        d_H *= ranfield_factor_hardening_modulus;

        sigma_y = ranfield_factor_yield_strength*GetYieldStrength(epsilon_p_eq2,d_sigma);
        d_sigma*= ranfield_factor_yield_strength;

        yield_condition = norm_dev - sqrt_2div3 * sigma_y;
    }

    if (i==100)
    {
        rLogger << "yield condition " << yield_condition << " delta_gamma " << delta_gamma << "\n";
        rLogger << "epsilon_p_eq " << rOldStaticData->mEpsilonPEq;
        rLogger << "total strain " << total_strain[0] << " " << total_strain[1] << " " << total_strain[2] << " " << total_strain[3] << " " << total_strain[4] << " " << total_strain[5] << " " <<"\n";
        rLogger << "plastic strain " << plastic_strain[0] << " " << plastic_strain[1] << " " << plastic_strain[2] << " " << plastic_strain[3] << " " << plastic_strain[4] << " " << plastic_strain[5] << " " <<"\n";
        rLogger << "back stress" << back_stress[0] << " " << back_stress[1] << " " << back_stress[2] << " " << back_stress[3] << " " << back_stress[4] << " " << back_stress[5] << " " <<"\n";
        throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::ReturnMapping3D] No convergence after 100 steps, check the source code.");
    }

    /* derivative of yield surface */
    df_dsigma[0]  = xi_trial[0]/norm_dev;
    df_dsigma[1]  = xi_trial[1]/norm_dev;
    df_dsigma[2]  = xi_trial[2]/norm_dev;
    df_dsigma[3]  = xi_trial[3]/norm_dev;
    df_dsigma[4]  = xi_trial[4]/norm_dev;
    df_dsigma[5]  = xi_trial[5]/norm_dev;

    //update static data
    if (rNewStaticData!=0)
    {
    	//update equivalent plastic strain
    	rNewStaticData->mEpsilonPEq = rOldStaticData->mEpsilonPEq + sqrt_2div3 * delta_gamma;

        //update backstress
        factor = sqrt_2div3 * (H2-H);
        rNewStaticData->mSigmaB[0] = rOldStaticData->mSigmaB[0] + factor*df_dsigma[0];
        rNewStaticData->mSigmaB[1] = rOldStaticData->mSigmaB[1] + factor*df_dsigma[1];
        rNewStaticData->mSigmaB[2] = rOldStaticData->mSigmaB[2] + factor*df_dsigma[2];
        rNewStaticData->mSigmaB[3] = rOldStaticData->mSigmaB[3] + factor*df_dsigma[3];
        rNewStaticData->mSigmaB[4] = rOldStaticData->mSigmaB[4] + factor*df_dsigma[4];
        rNewStaticData->mSigmaB[5] = rOldStaticData->mSigmaB[5] + factor*df_dsigma[5];

        //update plastic_strain
        rNewStaticData->mEpsilonP[0] = rOldStaticData->mEpsilonP[0] + delta_gamma*df_dsigma[0];
        rNewStaticData->mEpsilonP[1] = rOldStaticData->mEpsilonP[1] + delta_gamma*df_dsigma[1];
        rNewStaticData->mEpsilonP[2] = rOldStaticData->mEpsilonP[2] + delta_gamma*df_dsigma[2];
        rNewStaticData->mEpsilonP[3] = rOldStaticData->mEpsilonP[3] + 2.*delta_gamma*df_dsigma[3];  /* gamma */
        rNewStaticData->mEpsilonP[4] = rOldStaticData->mEpsilonP[4] + 2.*delta_gamma*df_dsigma[4];  /* gamma */
        rNewStaticData->mEpsilonP[5] = rOldStaticData->mEpsilonP[5] + 2.*delta_gamma*df_dsigma[5];  /* gamma */
    }

    //update stress
    if (rNewStress!=0)
    {
		factor  = 2.*mu*delta_gamma;
		factor2 = bulk_modulus*trace_epsilon;
		rNewStress->mEngineeringStress[0] = factor2+sigma_trial[0]-factor*df_dsigma[0];
		rNewStress->mEngineeringStress[1] = factor2+sigma_trial[1]-factor*df_dsigma[1];
		rNewStress->mEngineeringStress[2] = factor2+sigma_trial[2]-factor*df_dsigma[2];
		rNewStress->mEngineeringStress[3] =         sigma_trial[3]-factor*df_dsigma[3];
		rNewStress->mEngineeringStress[4] =         sigma_trial[4]-factor*df_dsigma[4];
		rNewStress->mEngineeringStress[5] =         sigma_trial[5]-factor*df_dsigma[5];
    }

    //update stiffness
    if (rNewTangent!=0)
    {
        double theta     = 1.-2.*mu*delta_gamma/norm_dev;
        double theta_bar = 1./(1.+(d_sigma+d_H)/(3.*mu))-(1.-theta);
        factor    = 2.*mu*theta;
        double factor_div3  = factor/3.;
        double factor_mul2div3  = 2.*factor_div3;
        double factor2   = -2.*mu*theta_bar;
        double factor3   = bulk_modulus;
        rNewTangent->mTangent[0]  =     (factor3 + factor_mul2div3 + factor2*df_dsigma[0]*df_dsigma[0]);
        rNewTangent->mTangent[1]  =     (factor3 - factor_div3     + factor2*df_dsigma[0]*df_dsigma[1]);
        rNewTangent->mTangent[2]  =     (factor3 - factor_div3     + factor2*df_dsigma[0]*df_dsigma[2]);
        rNewTangent->mTangent[3]  =     (							 factor2*df_dsigma[0]*df_dsigma[3]);
        rNewTangent->mTangent[4]  =     ( 						 factor2*df_dsigma[0]*df_dsigma[4]);
        rNewTangent->mTangent[5]  =     ( 						 factor2*df_dsigma[0]*df_dsigma[5]);
        rNewTangent->mTangent[6]  =     (rNewTangent->mTangent[1]);
        rNewTangent->mTangent[7]  =     (factor3 + factor_mul2div3 + factor2*df_dsigma[1]*df_dsigma[1]);
        rNewTangent->mTangent[8]  =     (factor3 - factor_div3     + factor2*df_dsigma[1]*df_dsigma[2]);
        rNewTangent->mTangent[9]  =     (						 factor2*df_dsigma[1]*df_dsigma[3]);
        rNewTangent->mTangent[10] =     (						 factor2*df_dsigma[1]*df_dsigma[4]);
        rNewTangent->mTangent[11] =     (						 factor2*df_dsigma[1]*df_dsigma[5]);
        rNewTangent->mTangent[12] =     (rNewTangent->mTangent[2]);
        rNewTangent->mTangent[13] =     (rNewTangent->mTangent[8]);
        rNewTangent->mTangent[14] =     (factor3 + factor_mul2div3 + factor2*df_dsigma[2]*df_dsigma[2]);
        rNewTangent->mTangent[15] =     (						 factor2*df_dsigma[2]*df_dsigma[3]);
        rNewTangent->mTangent[16] =     (						 factor2*df_dsigma[2]*df_dsigma[4]);
        rNewTangent->mTangent[17] =     (						 factor2*df_dsigma[2]*df_dsigma[5]);
        rNewTangent->mTangent[18] =     (rNewTangent->mTangent[3]);
        rNewTangent->mTangent[19] =     (rNewTangent->mTangent[9]);
        rNewTangent->mTangent[20] =     (rNewTangent->mTangent[15]);
        rNewTangent->mTangent[21] =     (	 0.5*factor 		 +factor2*df_dsigma[3]*df_dsigma[3]);
        rNewTangent->mTangent[22] =     (					      factor2*df_dsigma[3]*df_dsigma[4]);
        rNewTangent->mTangent[23] =     (					      factor2*df_dsigma[3]*df_dsigma[5]);
        rNewTangent->mTangent[24] =     (rNewTangent->mTangent[4]);
        rNewTangent->mTangent[25] =     (rNewTangent->mTangent[10]);
        rNewTangent->mTangent[26] =     (rNewTangent->mTangent[16]);
        rNewTangent->mTangent[27] =     (rNewTangent->mTangent[22]);
        rNewTangent->mTangent[28] =     (	 0.5*factor 		 +factor2*df_dsigma[4]*df_dsigma[4]);
        rNewTangent->mTangent[29] =     (					      factor2*df_dsigma[4]*df_dsigma[5]);
        rNewTangent->mTangent[30] =     (rNewTangent->mTangent[5]);
        rNewTangent->mTangent[31] =     (rNewTangent->mTangent[11]);
        rNewTangent->mTangent[32] =     (rNewTangent->mTangent[17]);
        rNewTangent->mTangent[33] =     (rNewTangent->mTangent[23]);
        rNewTangent->mTangent[34] =     (rNewTangent->mTangent[29]);
        rNewTangent->mTangent[35] =     (	 0.5*factor 		 +factor2*df_dsigma[5]*df_dsigma[5]);
    }
}


//! @brief ... calculates for a given equivalent plastic strain the radius of the yield surface
//! @param ... rEpsilonPEq equivalent plastic strain
//! @param ... rDSigmaDEpsilonP derivative of the yield strength with respect to the plastic strains (return value)
//! @return ... yield strength (radius of the yield surface)
double NuTo::ConstitutiveMisesPlasticity::GetYieldStrength(double rEpsilonPEq, double& rDSigmaDEpsilonP)const
{
    assert(mSigma.size()>0);
    std::vector<std::pair<double, double> >::const_iterator it(mSigma.begin());
    while (rEpsilonPEq>=it->first)
    {
        it++;
        if (it==mSigma.end())
        {
            // the maximum is reached, afterwards the yield strength remains constant
        	rDSigmaDEpsilonP = 0.;
            return (it-1)->second;
        }
    }

  	rDSigmaDEpsilonP = (it->second-(it-1)->second)/(it->first-(it-1)->first);
   	return rDSigmaDEpsilonP*(rEpsilonPEq-(it-1)->first)+(it-1)->second;
}

//! @brief ... calculates for a given equivalent plastic strain the hardening modulus
//! @param ... rEpsilonPEq equivalent plastic strain
//! @param ... rDSigmaDEpsilonP derivative of the hardening modulus with respect to the plastic strains (return value)
//! @return ... hardening modulus
double NuTo::ConstitutiveMisesPlasticity::GetHardeningModulus(double rEpsilonPEq, double& rDHDEpsilonP)const
{
    assert(mH.size()>0);
	std::vector<std::pair<double, double> >::const_iterator it(mH.begin());
    double H(0);
    if (mH.size()==1)
    {
        rDHDEpsilonP =it->second;
    	H=rEpsilonPEq*rDHDEpsilonP;
        return H;
    }
    else
    {
		do
		{
			it++;
			if (it==mH.end())
			{
				// the maximum is reached, afterwards use a constant slope
				it--;
				rDHDEpsilonP = it->second;
				H+=(rEpsilonPEq-it->first)*(it->second);
				return H;
			}
			if (rEpsilonPEq>it->first)
			{
				H+=((it)->first-(it-1)->first)*(it-1)->second;
			}
			else
			{
				it--;
				rDHDEpsilonP = it->second;
				H+=(rEpsilonPEq-it->first)*(it->second);
				return H;
			}
		}
		while(true);
    }
}

///////////////////////////////////////////////////////////////////////////

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::ConstitutiveMisesPlasticity::GetYoungsModulus() const
{
	return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::ConstitutiveMisesPlasticity::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::ConstitutiveMisesPlasticity::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::ConstitutiveMisesPlasticity::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}
//! @brief ... get initial yield strength
//! @return ... yield strength
double NuTo::ConstitutiveMisesPlasticity::GetInitialYieldStrength() const
{
	if(mSigma.size()==0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::GetInitialYieldStrength] Size of yield strength vector is zero.");
	return mSigma[0].second;
}

//! @brief ... set initial yield strength
//! @param rSigma ...  yield strength
void NuTo::ConstitutiveMisesPlasticity::SetInitialYieldStrength(double rSigma)
{
	if(mSigma.size()==0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::SetInitialYieldStrength] Size of yield strength vector is zero.");
	if (rSigma<=0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::SetInitialYieldStrength] Initial yield strength has to be positive.");
	mSigma[0].second = rSigma;
    this->SetParametersValid();
}

//! @brief ... get yield strength for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding yield strength
NuTo::FullMatrix<double> NuTo::ConstitutiveMisesPlasticity::GetYieldStrength() const
{
	NuTo::FullMatrix<double> returnMatrix(mSigma.size(),2);
	for (unsigned int count=0; count<mSigma.size(); count++)
	{
		returnMatrix(count,0) = mSigma[count].first;
		returnMatrix(count,1) = mSigma[count].second;
	}
	return returnMatrix;
}

//! @brief ... add yield strength
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  yield strength
void NuTo::ConstitutiveMisesPlasticity::AddYieldStrength(double rEpsilon, double rSigma)
{
	if (rEpsilon<=0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AddYieldStrength] Equivalente strain has to be positive.");
	if (rSigma<=0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AddYieldStrength] Yield strength has to be positive.");
	std::vector<std::pair<double, double> >::iterator it;
	for (it = mSigma.begin(); it!=mSigma.end(); it++)
	{
		if (it->first>rEpsilon)
		{
			if (it!=mSigma.begin())
			{
				if ((it-1)->second>(it->second))
					throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AddYieldStrength] The yield strength can only increase for increasing epsilon equivalent.");
			}
			break;
		}
	}
	mSigma.insert(it,1,std::pair<double,double>(rEpsilon,rSigma));
    this->SetParametersValid();
}

//! @brief ... get initial hardening modulus
//! @return ... hardening modulus
double NuTo::ConstitutiveMisesPlasticity::GetInitialHardeningModulus() const
{
	if(mH.size()==0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::SetInitialHardeningModulus] Size of hardening modulus vector is zero.");
	return mH[0].second;
}

//! @brief ... set initial hardening modulus
//! @param rH ...  hardening modulus
void NuTo::ConstitutiveMisesPlasticity::SetInitialHardeningModulus(double rH)
{
	if(mH.size()==0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::SetInitialHardeningModulus] Size of hardening modulus vector is zero.");
	if (rH<0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::SetInitialHardeningModulus] Initial hardening modulus must not be negative.");
	mH[0].second = rH;
    this->SetParametersValid();
}

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double> NuTo::ConstitutiveMisesPlasticity::GetHardeningModulus() const
{
	NuTo::FullMatrix<double> returnMatrix(mH.size(),2);
	for (unsigned int count=0; count<mH.size(); count++)
	{
		returnMatrix(count,0) = mH[count].first;
		returnMatrix(count,1) = mH[count].second;
	}
	return returnMatrix;
}

//! @brief ... add hardening modulus
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  hardening modulus
void NuTo::ConstitutiveMisesPlasticity::AddHardeningModulus(double rEpsilon, double rH)
{
	if (rEpsilon<=0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AddHardeningModulus] Equivalente strain has to be positive.");
	if (rH<0)
		throw MechanicsException("[NuTo::ConstitutiveMisesPlasticity::AddHardeningModulus] Hardening modul must not be negative.");
	std::vector<std::pair<double, double> >::iterator it;
	for (it = mH.begin(); it!=mH.end(); it++)
	{
		if (it->first>rEpsilon)
		{
			break;
		}
	}
	mH.insert(it,1,std::pair<double,double>(rEpsilon,rH));
    this->SetParametersValid();
}

///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::ConstitutiveMisesPlasticity::GetType() const
{
    return NuTo::Constitutive::LINEAR_ELASTIC;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::ConstitutiveMisesPlasticity::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::TETRAHEDRON4N:
        return true;
    case NuTo::Element::TETRAHEDRON10N:
        return true;
    case NuTo::Element::BRICK8N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::ConstitutiveMisesPlasticity::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::ConstitutiveMisesPlasticity::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check yield strength is positive
//! @param rSigma ... yield strength
void NuTo::ConstitutiveMisesPlasticity::CheckYieldStrength(std::vector<std::pair<double, double> > rSigma) const
{
	if (rSigma.size()==0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckYieldStrength] At least an initial yield strength is required.");
    }
	if (rSigma[0].second <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckYieldStrength] The initial yield strength must be a positive value.");
    }
	rSigma[0].first = 0.;

	for (unsigned int count=1; count<rSigma.size(); count++)
	{
		if (rSigma[count-1].first>=rSigma[count].first)
	        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckYieldStrength] For multilinear plasticity, the epsilon should always increase.");
		if (rSigma[count-1].second>rSigma[count].second)
	        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckYieldStrength] For multilinear plasticity, the yield strength should always increase.");
	}
}

//! @brief ... check hardening modulus
//! @param rH ... hardening modulus
void NuTo::ConstitutiveMisesPlasticity::CheckHardeningModulus(std::vector<std::pair<double, double> > rH) const
{
	if (rH.size()==0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckHardeningModulus] At least an initial hardening modulus is required.");
    }
	if (rH[0].second < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckHardeningModulus] The initial hardening modulus must not be a negative value.");
    }
	rH[0].first = 0.;

	for (unsigned int count=1; count<rH.size(); count++)
	{
		if (rH[count-1].first>=rH[count].first)
	        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckHardeningModulus] For multilinear plasticity, the epsilon should always increase.");
		if (rH[count].second<0)
	        throw NuTo::MechanicsException("[NuTo::ConstitutiveMisesPlasticity::CheckHardeningModulus] For multilinear plasticity, the hardening modulus should always be positive.");
	}
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstitutiveMisesPlasticity::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus: " << this->mE << "\n";
    rLogger << "    Poisson's ratio: " << this->mNu << "\n";
	rLogger << "    multilinear yield strength: (interval,epsilon,sigma)" << "\n";
    for (unsigned int count=0; count<mSigma.size(); count++)
    	rLogger << "       " << count<< " : " << this->mSigma[count].first << "    " << this->mSigma[count].second << "\n";
	rLogger << "    multilinear hardening modulus: (interval,epsilon,H')" << "\n";
    for (unsigned int count=0; count<mH.size(); count++)
    	rLogger << "       " << count<< " : " << this->mH[count].first << "    " << this->mH[count].second << "\n";
}

// check parameters
void NuTo::ConstitutiveMisesPlasticity::CheckParameters()const
{
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckYieldStrength(this->mSigma);
    this->CheckHardeningModulus(this->mH);
}
