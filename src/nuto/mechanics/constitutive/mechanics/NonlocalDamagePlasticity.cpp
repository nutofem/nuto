// $ld: $ 
// NonlocalDamagePlasticity.cpp
// created Apr 26, 2010 by Joerg F. Unger
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticity.h"
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
#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

NuTo::NonlocalDamagePlasticity::NonlocalDamagePlasticity() : ConstitutiveEngineeringStressStrain()
{
	std::cout << "[NuTo::NonlocalDamagePlasticity::NonlocalDamagePlasticity]" << std::endl;
	mE = 0.;
	mNu = 0.;
	mNonlocalRadius = 1.;
	mTensileStrength = 0.;
	mCompressiveStrength = 0.;
	mBiaxialCompressiveStrength = 0.;
	mFractureEnergy = 0.;
	SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NonlocalDamagePlasticity::serialize(Archive & ar, const unsigned int version)
    {
        std::cout << "start serialization of linear elastic" << std::endl;
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
           & BOOST_SERIALIZATION_NVP(mE)
           & BOOST_SERIALIZATION_NVP(mNu)
           & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
           & BOOST_SERIALIZATION_NVP(mTensileStrength)
           & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
           & BOOST_SERIALIZATION_NVP(mBiaxialCompressiveStrength)
           & BOOST_SERIALIZATION_NVP(mFractureEnergy);
        std::cout << "finish serialization of linear elastic" << std::endl;
    }
#endif // ENABLE_SERIALIZATION

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::NonlocalDamagePlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	//TODO check this routine if it is necessary
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::NonlocalDamagePlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	//TODO check this routine if it is necessary
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::NonlocalDamagePlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    //TODO check this routine if it is necessary
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Material model not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... engineering stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
	      const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Material model not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    assert(rElement->GetSection()!=0);
    if (rElement->GetSection()->GetType()==Section::PLANE_STRAIN)
    {
    	ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain *oldStaticData =
    			dynamic_cast<ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain*>(rElement->GetStaticData(rIp));
        if (oldStaticData==0)
        	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Static data of element is not of type ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.");
    	// no new history variables, since no update has to be performed
        this->ReturnMapping2D(rElement,engineeringStrain,oldStaticData, 0);

		// calculate Engineering stress
		//rEngineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * engineeringStrain.mEngineeringStrain[1];
		//rEngineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * engineeringStrain.mEngineeringStrain[0];
		//rEngineeringStress.mEngineeringStress[2] = C33 * engineeringStrain.mEngineeringStrain[2] ;

    }
    else
    {
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Plane stress is not implemented.");
    }
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
	// check if parameters are valid
	if (this->mParametersValid == false)
	{
		//throw an exception giving information related to the wrong parameter
	    CheckParameters();
	    //if there is no exception thrown there is a problem with the source code
	    //since every time a material parameter is changed, the parametes should be checked
	    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
	}
	// calculate engineering strain
	EngineeringStrain2D engineeringStrain;
	rDeformationGradient.GetEngineeringStrain(engineeringStrain);

	assert(rElement->GetSection()!=0);
	if (rElement->GetSection()->GetType()==Section::PLANE_STRAIN)
	{
	    // calculate coefficients of the material matrix
	    double C11, C12, C33;
	    this->CalculateCoefficients3D(C11, C12, C33);

	    // calculate Engineering stress
	    rEngineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * engineeringStrain.mEngineeringStrain[1];
	    rEngineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * engineeringStrain.mEngineeringStrain[0];
	    rEngineeringStress.mEngineeringStress[2] = C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
	    rEngineeringStress.mEngineeringStress[3] = C33 * engineeringStrain.mEngineeringStrain[2] ;
	    rEngineeringStress.mEngineeringStress[4] = 0.;
	    rEngineeringStress.mEngineeringStress[5] = 0.;
	}
	else
	{
	    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
	}
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain3D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    rEngineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * (engineeringStrain.mEngineeringStrain[1]+engineeringStrain.mEngineeringStrain[2]);
    rEngineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[2]);
    rEngineeringStress.mEngineeringStress[2] = C11 * engineeringStrain.mEngineeringStrain[2] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
    rEngineeringStress.mEngineeringStress[3] = C44 * engineeringStrain.mEngineeringStrain[3] ;
    rEngineeringStress.mEngineeringStress[4] = C44 * engineeringStrain.mEngineeringStrain[4] ;
    rEngineeringStress.mEngineeringStress[5] = C44 * engineeringStrain.mEngineeringStrain[5] ;
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient,
		ConstitutiveTangentLocal1x1& rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // store tangent at the output object
    rTangent.mTangent = mE;
    rTangent.SetSymmetry(true);
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient,
		ConstitutiveTangentLocal3x3& rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

	const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
	if (theSection->GetType()==Section::PLANE_STRAIN)
	{
	    // calculate coefficients of the material matrix
	    double C11, C12, C33;
	    this->CalculateCoefficients3D(C11, C12, C33);

	    // store tangent at the output object
	     rTangent.mTangent[ 0] = C11;
	     rTangent.mTangent[ 1] = C12;
	     rTangent.mTangent[ 2] = 0;

	     rTangent.mTangent[ 3] = C12;
	     rTangent.mTangent[ 4] = C11;
	     rTangent.mTangent[ 5] = 0;

	     rTangent.mTangent[ 6] = 0.;
	     rTangent.mTangent[ 7] = 0.;
	     rTangent.mTangent[ 8] = C33;
	}
	else
	{
	    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
	}

    rTangent.SetSymmetry(true);
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient,
		ConstitutiveTangentLocal6x6& rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // store tangent at the output object
    rTangent.mTangent[ 0] = C11;
    rTangent.mTangent[ 1] = C12;
    rTangent.mTangent[ 2] = C12;
    rTangent.mTangent[ 3] = 0;
    rTangent.mTangent[ 4] = 0;
    rTangent.mTangent[ 5] = 0;

    rTangent.mTangent[ 6] = C12;
    rTangent.mTangent[ 7] = C11;
    rTangent.mTangent[ 8] = C12;
    rTangent.mTangent[ 9] = 0;
    rTangent.mTangent[10] = 0;
    rTangent.mTangent[11] = 0;

    rTangent.mTangent[12] = C12;
    rTangent.mTangent[13] = C12;
    rTangent.mTangent[14] = C11;
    rTangent.mTangent[15] = 0;
    rTangent.mTangent[16] = 0;
    rTangent.mTangent[17] = 0;

    rTangent.mTangent[18] = 0;
    rTangent.mTangent[19] = 0;
    rTangent.mTangent[20] = 0;
    rTangent.mTangent[21] = C44;
    rTangent.mTangent[22] = 0;
    rTangent.mTangent[23] = 0;

    rTangent.mTangent[24] = 0;
    rTangent.mTangent[25] = 0;
    rTangent.mTangent[26] = 0;
    rTangent.mTangent[27] = 0;
    rTangent.mTangent[28] = C44;
    rTangent.mTangent[29] = 0;

    rTangent.mTangent[30] = 0;
    rTangent.mTangent[31] = 0;
    rTangent.mTangent[32] = 0;
    rTangent.mTangent[33] = 0;
    rTangent.mTangent[34] = 0;
    rTangent.mTangent[35] = C44;

    rTangent.SetSymmetry(true);
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::NonlocalDamagePlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::NonlocalDamagePlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::NonlocalDamagePlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
	   //no static data required -> empty routine
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
		const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for 1D.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
		const ElementBase* rElement) const
{
	if (rElement->GetSection()->GetType()==NuTo::Section::PLANE_STRESS)
		throw MechanicsException("[NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for plane stress.");
	else
	    return new ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain();
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
		const ElementBase* rElement) const
{
	return new ConstitutiveStaticDataNonlocalDamagePlasticity3D();
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate engineering strain
    EngineeringStrain1D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);
    return 0.5 * engineeringStrain.mEngineeringStrain * this->mE * engineeringStrain.mEngineeringStrain;
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    EngineeringStress2D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

	const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
	if (theSection->GetType()==Section::PLANE_STRAIN)
	{
		// calculate coefficients of the material matrix
		double C11, C12, C33;
		this->CalculateCoefficients3D(C11, C12, C33);

		// calculate Engineering stress
		engineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * engineeringStrain.mEngineeringStrain[1];
		engineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * engineeringStrain.mEngineeringStrain[0];
		engineeringStress.mEngineeringStress[2] = C33 * engineeringStrain.mEngineeringStrain[2] ;
	}
	else
	{
		throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Plane stress is to be implemented.");
	}
    return 0.5*(
    		engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]);
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain3D engineeringStrain;
    EngineeringStress3D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    engineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * (engineeringStrain.mEngineeringStrain[1]+engineeringStrain.mEngineeringStrain[2]);
    engineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[2]);
    engineeringStress.mEngineeringStress[2] = C11 * engineeringStrain.mEngineeringStrain[2] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
    engineeringStress.mEngineeringStress[3] = C44 * engineeringStrain.mEngineeringStrain[3] ;
    engineeringStress.mEngineeringStress[4] = C44 * engineeringStrain.mEngineeringStrain[4] ;
    engineeringStress.mEngineeringStress[5] = C44 * engineeringStrain.mEngineeringStrain[5] ;

    return 0.5*(
    		engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]
		   +engineeringStrain.mEngineeringStrain[3]*engineeringStress.mEngineeringStress[3]
		   +engineeringStrain.mEngineeringStrain[4]*engineeringStress.mEngineeringStress[4]
		   +engineeringStrain.mEngineeringStrain[5]*engineeringStress.mEngineeringStress[5]);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
	return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
	return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
	return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

// calculate coefficients of the material matrix
void NuTo::NonlocalDamagePlasticity::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2*(1.0 + this->mNu));
}

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::NonlocalDamagePlasticity::GetYoungsModulus() const
{
	return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::NonlocalDamagePlasticity::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::NonlocalDamagePlasticity::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::NonlocalDamagePlasticity::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

//! @brief ... get nonlocal radius
//! @return ... nonlocal radius
double NuTo::NonlocalDamagePlasticity::GetNonlocalRadius() const
{
    return mNonlocalRadius;
}

//! @brief ... set nonlocal radius
//! @param rRadius...  nonlocal radius
void NuTo::NonlocalDamagePlasticity::SetNonlocalRadius(double rNonlocalRadius)
{
    this->CheckNonlocalRadius(rNonlocalRadius);
    this->mNonlocalRadius = rNonlocalRadius;
    this->SetParametersValid();
}
//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::NonlocalDamagePlasticity::GetTensileStrength() const
{
    return mTensileStrength;
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::NonlocalDamagePlasticity::SetTensileStrength(double rTensileStrength)
{
    this->CheckNonlocalRadius(rTensileStrength);
    this->mTensileStrength = rTensileStrength;
    this->SetParametersValid();
}

//! @brief ... get compressive strength
//! @return ... compressive strength
double NuTo::NonlocalDamagePlasticity::GetCompressiveStrength() const
{
    return mCompressiveStrength;
}

//! @brief ... set compressive strength
//! @param rCompressiveStrength...  compressive strength
void NuTo::NonlocalDamagePlasticity::SetCompressiveStrength(double rCompressiveStrength)
{
    this->CheckNonlocalRadius(rCompressiveStrength);
    this->mCompressiveStrength = rCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get biaxial compressive strength
//! @return ... biaxial compressive strength
double NuTo::NonlocalDamagePlasticity::GetBiaxialCompressiveStrength() const
{
    return mBiaxialCompressiveStrength;
}

//! @brief ... set biaxial compressive strength
//! @param rBiaxialCompressiveStrength...  biaxial compressive strength
void NuTo::NonlocalDamagePlasticity::SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength)
{
    this->CheckNonlocalRadius(rBiaxialCompressiveStrength);
    this->mBiaxialCompressiveStrength = rBiaxialCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get fracture energy
//! @return ... fracture energy
double NuTo::NonlocalDamagePlasticity::GetFractureEnergy() const
{
    return mFractureEnergy;
}

//! @brief ... set fracture energy
//! @param rFractureEnergy... fracture energy
void NuTo::NonlocalDamagePlasticity::SetFractureEnergy(double rFractureEnergy)
{
    this->CheckNonlocalRadius(rFractureEnergy);
    this->mFractureEnergy = rFractureEnergy;
    this->SetParametersValid();
}
///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::NonlocalDamagePlasticity::GetType() const
{
    return NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::NonlocalDamagePlasticity::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return true;
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
        return true;
    case NuTo::Element::PLANE2D6N:
        return true;
    case NuTo::Element::TETRAHEDRON4N:
        return true;
    case NuTo::Element::TETRAHEDRON10N:
        return true;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::NonlocalDamagePlasticity::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::NonlocalDamagePlasticity::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::NonlocalDamagePlasticity::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}
//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::NonlocalDamagePlasticity::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::NonlocalDamagePlasticity::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if biaxial compressive strength is positive
//! @param rBiaxialCompressiveStrength ... biaxial compressive strength
void NuTo::NonlocalDamagePlasticity::CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const
{
    if (rBiaxialCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be a positive value.");
    }
    if (rBiaxialCompressiveStrength <= mCompressiveStrength)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be higher than the uniaxial compressive strength.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::NonlocalDamagePlasticity::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}
//! @brief ... performs the return mapping procedure for the plasticity model
void NuTo::NonlocalDamagePlasticity::ReturnMapping2D(
		const ElementBase* rElement,
		const EngineeringStrain2D& rStrain,
		const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* rOldStaticData,
		ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* rNewStaticData)const
{
    //TODO
	double minCutbackFactor(1e-3);
	double stress[4];
    double mM(1);

    double l_element(sqrt(rElement->CalculateArea()*0.5));

    double kappa_d = 15.*mFractureEnergy/(16.*mTensileStrength);
    double l_eq_plane = mNonlocalRadius;
    double l_eq_circ  = mNonlocalRadius*rElement->GetSection()->GetThickness()/l_element;

    double e_mod = mE; //modify that one in the case of random fields
    double   nu  = mNu;
    double f_ct  = mTensileStrength;
    double f_ck  = mCompressiveStrength;
    double f_ck2 = mBiaxialCompressiveStrength;
    double     G = mFractureEnergy;

    assert(f_ck2>f_ck);


	double curTotalStrain[4];
	double curPlastStrain[4];
	double lastPlastStrain[4];
	double deltaStrain[4];

	/* apply strains in steps */
	deltaStrain[0] = rStrain.mEngineeringStrain[0]-rOldStaticData->GetPrevStrain()[0];
	deltaStrain[1] = rStrain.mEngineeringStrain[0]-rOldStaticData->GetPrevStrain()[1];
	deltaStrain[2] = rStrain.mEngineeringStrain[0]-rOldStaticData->GetPrevStrain()[2];
	deltaStrain[3] = 0;

    lastPlastStrain[0] =  rOldStaticData->mEpsilonP[0];
    lastPlastStrain[1] =  rOldStaticData->mEpsilonP[1];
    lastPlastStrain[2] =  rOldStaticData->mEpsilonP[2];
    lastPlastStrain[3] =  rOldStaticData->mEpsilonP[3];

    curPlastStrain[0] =  lastPlastStrain[0];
    curPlastStrain[1] =  lastPlastStrain[1];
    curPlastStrain[2] =  lastPlastStrain[2];
    curPlastStrain[3] =  lastPlastStrain[3];

	double cutbackFactor = 1, deltaCutbackFactor=1.;
	bool finished(false);
	bool converged;
	while (cutbackFactor>minCutbackFactor && !finished)
	{

		curTotalStrain[0] = rOldStaticData->mPrevStrain[0]+cutbackFactor*deltaStrain[0];
		curTotalStrain[1] = rOldStaticData->mPrevStrain[1]+cutbackFactor*deltaStrain[1];
		curTotalStrain[2] = rOldStaticData->mPrevStrain[2]+cutbackFactor*deltaStrain[2];
		curTotalStrain[3] = 0.;

        try
        {
        	converged = false;
        	// perform return mapping

        }
        catch (NuTo::MechanicsException& e)
        {
        	e.AddMessage("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] Error performing return mapping procedure.");
        	throw e;
        }
        catch (...)
        {
           	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] Error performing return mapping procedure.");
        }

        if (converged)
        {
        	if (cutbackFactor==1.)
        		finished=true;
        	else
        	{
        	    curPlastStrain[0] =  lastPlastStrain[0];
        	    curPlastStrain[1] =  lastPlastStrain[1];
        	    curPlastStrain[2] =  lastPlastStrain[2];
        	    curPlastStrain[3] =  lastPlastStrain[3];
        	    cutbackFactor +=cutbackFactor;
				if (cutbackFactor>1)
				{
					deltaCutbackFactor -= cutbackFactor -1.;
					cutbackFactor = 1.;
				}
        	}
        }
        else
        {
			/* no convergence in return mapping*/
        	deltaCutbackFactor*=0.5;
			deltaCutbackFactor-=deltaCutbackFactor;
        }
	}
	if (deltaCutbackFactor<=minCutbackFactor)
	{
		throw MechanicsException("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] No convergence can be obtained in the return mapping procedure.");
	}

	//here set new stress, tangent and static data
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::NonlocalDamagePlasticity::Info(unsigned short rVerboseLevel) const
{
    this->ConstitutiveBase::Info(rVerboseLevel);
    std::cout << "    Young's modulus      : " << this->mE << std::endl;
    std::cout << "    Poisson's ratio      : " << this->mNu << std::endl;
    std::cout << "    nonlocal radius      : " << this->mNonlocalRadius << std::endl;
    std::cout << "    tensile strength     : " << this->mTensileStrength << std::endl;
    std::cout << "    compressive strength : " << this->mCompressiveStrength << std::endl;
    std::cout << "    biaxial compressive strength : " << this->mBiaxialCompressiveStrength << std::endl;
    std::cout << "    fracture energy      : " << this->mFractureEnergy << std::endl;
}

// check parameters
void NuTo::NonlocalDamagePlasticity::CheckParameters()const
{
    this->CheckBiaxialCompressiveStrength(this->mNu);
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckNonlocalRadius(this->mNu);
    this->CheckTensileStrength(this->mNu);
    this->CheckCompressiveStrength(this->mNu);
    this->CheckFractureEnergy(this->mNu);
}


