// $Id: LinearElastic1D.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/LinearElastic.h"
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
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

NuTo::LinearElastic::LinearElastic() : ConstitutiveEngineeringStressStrain(), ConstitutivePiolaKirchhoffIIGreenLagrange()
{
	mE = 0.;
	mNu = 0.;
	mRho = 0.;
	SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::LinearElastic::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
   std::cout << "start serialize LinearElastic" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
      & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutivePiolaKirchhoffIIGreenLagrange)
      & BOOST_SERIALIZATION_NVP(mE)
      & BOOST_SERIALIZATION_NVP(mNu);
#ifdef DEBUG_SERIALIZATION
   std::cout << "finish serialize LinearElastic" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LinearElastic)
#endif // ENABLE_SERIALIZATION

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::LinearElastic::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	rEngineeringPlasticStrain.mEngineeringStrain[0] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[1] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[2] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[3] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[4] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[5] = 0.;
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::LinearElastic::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	rEngineeringPlasticStrain.mEngineeringStrain[0] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[1] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[2] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[3] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[4] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[5] = 0.;
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::LinearElastic::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
								  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
	rEngineeringPlasticStrain.mEngineeringStrain[0] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[1] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[2] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[3] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[4] = 0.;
	rEngineeringPlasticStrain.mEngineeringStrain[5] = 0.;
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain1D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // calculate Engineering stress
    rEngineeringStress.mEngineeringStress = mE * engineeringStrain.mEngineeringStrain;
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... engineering stress
void NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
	      const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
	// check if parameters are valid
	if (this->mParametersValid == false)
	{
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
	    throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
	}
	// calculate engineering strain
	EngineeringStrain1D engineeringStrain;
	rDeformationGradient.GetEngineeringStrain(engineeringStrain);

	// calculate Engineering stress
	rEngineeringStress.mEngineeringStress[0] = mE * engineeringStrain.mEngineeringStrain;
	rEngineeringStress.mEngineeringStress[1] = 0.;
	rEngineeringStress.mEngineeringStress[2] = 0.;
	rEngineeringStress.mEngineeringStress[3] = 0.;
	rEngineeringStress.mEngineeringStress[4] = 0.;
	rEngineeringStress.mEngineeringStress[5] = 0.;
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
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
		rEngineeringStress.mEngineeringStress[2] = C33 * engineeringStrain.mEngineeringStrain[2] ;
    }
    else
    {
    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
    }
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
	// check if parameters are valid
	if (this->mParametersValid == false)
	{
		//throw an exception giving information related to the wrong parameter
	    CheckParameters();
	    //if there is no exception thrown there is a problem with the source code
	    //since every time a material parameter is changed, the parametes should be checked
	    throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
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
	    throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
	}
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
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
//  Damage /////////////////////////////////////
 //! @brief ... calculate isotropic damage from deformation gradient in 1D
 //! @param rElement ... element
 //! @param rIp ... integration point
 //! @param rDeformationGradient ... deformation gradient
 //! @param rDamage ... damage variable
 void NuTo::LinearElastic::GetDamage(const ElementBase* rElement, int rIp,
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
 void NuTo::LinearElastic::GetDamage(const ElementBase* rElement, int rIp,
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
 void NuTo::LinearElastic::GetDamage(const ElementBase* rElement, int rIp,
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
void NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient,
		ConstitutiveTangentBase* rTangent) const
{
	ConstitutiveTangentLocal1x1 *tangent(rTangent->AsConstitutiveTangentLocal1x1());
	// check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // store tangent at the output object
    tangent->mTangent = mE;
    tangent->SetSymmetry(true);
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient,
		ConstitutiveTangentBase* rTangent) const
{
	ConstitutiveTangentLocal3x3 *tangent(rTangent->AsConstitutiveTangentLocal3x3());
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

	const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
	if (theSection->GetType()==Section::PLANE_STRAIN)
	{
	    // calculate coefficients of the material matrix
	    double C11, C12, C33;
	    this->CalculateCoefficients3D(C11, C12, C33);

	    // store tangent at the output object
	    tangent->mTangent[ 0] = C11;
	    tangent->mTangent[ 1] = C12;
	    tangent->mTangent[ 2] = 0;

	    tangent->mTangent[ 3] = C12;
	    tangent->mTangent[ 4] = C11;
	    tangent->mTangent[ 5] = 0;

	    tangent->mTangent[ 6] = 0.;
	    tangent->mTangent[ 7] = 0.;
	    tangent->mTangent[ 8] = C33;
	}
	else
	{
	    throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
	}

    rTangent->SetSymmetry(true);
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient,
		ConstitutiveTangentBase* rTangent) const
{
	ConstitutiveTangentLocal6x6 *tangent(rTangent->AsConstitutiveTangentLocal6x6());

	// check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    double *data(tangent->mTangent);

    // store tangent at the output object
    data[ 0] = C11;
    data[ 1] = C12;
    data[ 2] = C12;
    data[ 3] = 0;
    data[ 4] = 0;
    data[ 5] = 0;

    data[ 6] = C12;
    data[ 7] = C11;
    data[ 8] = C12;
    data[ 9] = 0;
    data[10] = 0;
    data[11] = 0;

    data[12] = C12;
    data[13] = C12;
    data[14] = C11;
    data[15] = 0;
    data[16] = 0;
    data[17] = 0;

    data[18] = 0;
    data[19] = 0;
    data[20] = 0;
    data[21] = C44;
    data[22] = 0;
    data[23] = 0;

    data[24] = 0;
    data[25] = 0;
    data[26] = 0;
    data[27] = 0;
    data[28] = C44;
    data[29] = 0;

    data[30] = 0;
    data[31] = 0;
    data[32] = 0;
    data[33] = 0;
    data[34] = 0;
    data[35] = C44;

    tangent->SetSymmetry(true);
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
	   //no static data required -> empty routine
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
	//no need to update tmp static data
}


//! @brief ... update mp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
	//no need to update tmp static data
}

//! @brief ... update mp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
	//no need to update tmp static data
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::LinearElastic::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
    return 0;
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::LinearElastic::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
	return 0;
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::LinearElastic::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
	return 0;
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
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
double NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    EngineeringStress2D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

	const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
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
		throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain] Plane stress is to be implemented.");
	}
/*#ifdef DEBUG
    std::cout << "strain ";
    std::cout << engineeringStrain.mEngineeringStrain[0] << " ";
    std::cout << engineeringStrain.mEngineeringStrain[1] << " ";
    std::cout << engineeringStrain.mEngineeringStrain[2]<< std::endl;
    std::cout << "stress ";
    std::cout << engineeringStress.mEngineeringStress[0] << " ";
    std::cout << engineeringStress.mEngineeringStress[1] << " ";
    std::cout << engineeringStress.mEngineeringStress[2] << std::endl;
    std::cout << "energy " << 0.5*(
            engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]) << std::endl<<std::endl;
#endif //DEBUG
*/
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
double NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
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
double NuTo::LinearElastic::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
	return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
	return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
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
void NuTo::LinearElastic::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::LinearElastic::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::LinearElastic::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

///////////////////////////////////////////////////////////////////////////

// Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
//! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
void NuTo::LinearElastic::GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient1D& rDeformationGradient, SecondPiolaKirchhoffStress1D& rSecondPiolaKirchhoffStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain] Check the material parameters.");
    }

    // calculate engineering strain
    GreenLagrangeStrain1D greenLagrangeStrain;
    rDeformationGradient.GetGreenLagrangeStrain(greenLagrangeStrain);

    // calculate second piola kirchhoff stress
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress = mE * greenLagrangeStrain.mGreenLagrangeStrain;

}


// Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
//! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
void NuTo::LinearElastic::GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient2D& rDeformationGradient, SecondPiolaKirchhoffStress2D& rSecondPiolaKirchhoffStress) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain] To be implemented.");
}


// Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
//! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
void NuTo::LinearElastic::GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
			  const DeformationGradient3D& rDeformationGradient, SecondPiolaKirchhoffStress3D& rSecondPiolaKirchhoffStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain] Check the material parameters.");
    }
    // calculate engineering strain
    GreenLagrangeStrain3D greenLagrangeStrain;
    rDeformationGradient.GetGreenLagrangeStrain(greenLagrangeStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[0] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[0] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[1]+greenLagrangeStrain.mGreenLagrangeStrain[2]);
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[1] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[1] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[0]+greenLagrangeStrain.mGreenLagrangeStrain[2]);
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[2] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[2] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[0]+greenLagrangeStrain.mGreenLagrangeStrain[1]);
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[3] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[3] ;
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[4] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[4] ;
    rSecondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[5] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[5] ;

}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::LinearElastic::GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
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
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain] Check the material parameters.");
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
void NuTo::LinearElastic::GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient,
		ConstitutiveTangentLocal3x3& rTangent) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain] To be implemented.");
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::LinearElastic::GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
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
    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain] Check the material parameters.");
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
void NuTo::LinearElastic::UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain(ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain(ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
	//no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::LinearElastic::UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain(ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
	//no static data required -> empty routine
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::LinearElastic::AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain1D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
    return 0;
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::LinearElastic::AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain2D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
    return 0;

}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::LinearElastic::AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain3D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
    return 0;
}



//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetTotalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain] Check the material parameters.");
    }

    // calculate engineering strain
    GreenLagrangeStrain1D greenLagrangeStrain;
    rDeformationGradient.GetGreenLagrangeStrain(greenLagrangeStrain);
    return 0.5 * greenLagrangeStrain.mGreenLagrangeStrain * this->mE * greenLagrangeStrain.mGreenLagrangeStrain;
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetTotalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain] To be implemented.");
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetTotalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetTotalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain] Check the material parameters.");
    }
    // calculate engineering strain
    GreenLagrangeStrain3D greenLagrangeStrain;
    SecondPiolaKirchhoffStress3D secondPiolaKirchhoffStress;
    rDeformationGradient.GetGreenLagrangeStrain(greenLagrangeStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[0] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[0] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[1]+greenLagrangeStrain.mGreenLagrangeStrain[2]);
    secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[1] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[1] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[0]+greenLagrangeStrain.mGreenLagrangeStrain[2]);
    secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[2] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[2] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[0]+greenLagrangeStrain.mGreenLagrangeStrain[1]);
    secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[3] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[3] ;
    secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[4] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[4] ;
    secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[5] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[5] ;

    return 0.5*(
    		greenLagrangeStrain.mGreenLagrangeStrain[0]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[0]
           +greenLagrangeStrain.mGreenLagrangeStrain[1]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[1]
           +greenLagrangeStrain.mGreenLagrangeStrain[2]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[2]
		   +greenLagrangeStrain.mGreenLagrangeStrain[3]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[3]
		   +greenLagrangeStrain.mGreenLagrangeStrain[4]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[4]
		   +greenLagrangeStrain.mGreenLagrangeStrain[5]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[5]);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient1D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain] Check the material parameters.");
    }

    // calculate engineering strain
    GreenLagrangeStrain1D greenLagrangeStrain;
    rDeformationGradient.GetGreenLagrangeStrain(greenLagrangeStrain);
    return 0.5 * greenLagrangeStrain.mGreenLagrangeStrain * this->mE * greenLagrangeStrain.mGreenLagrangeStrain;
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient2D& rDeformationGradient) const
{
	throw MechanicsException("[NuTo::LinearElastic::GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain] To be implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::LinearElastic::GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
		const DeformationGradient3D& rDeformationGradient) const
{
	// check if parameters are valid
	if (this->mParametersValid == false)
	{
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::LinearElastic::GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain] Check the material parameters.");
	}
	// calculate engineering strain
	GreenLagrangeStrain3D greenLagrangeStrain;
	SecondPiolaKirchhoffStress3D secondPiolaKirchhoffStress;
	rDeformationGradient.GetGreenLagrangeStrain(greenLagrangeStrain);

	// calculate coefficients of the material matrix
	double C11, C12, C44;
	this->CalculateCoefficients3D(C11, C12, C44);

	// calculate Engineering stress
	secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[0] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[0] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[1]+greenLagrangeStrain.mGreenLagrangeStrain[2]);
	secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[1] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[1] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[0]+greenLagrangeStrain.mGreenLagrangeStrain[2]);
	secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[2] = C11 * greenLagrangeStrain.mGreenLagrangeStrain[2] + C12 * (greenLagrangeStrain.mGreenLagrangeStrain[0]+greenLagrangeStrain.mGreenLagrangeStrain[1]);
	secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[3] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[3] ;
	secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[4] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[4] ;
	secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[5] = C44 * greenLagrangeStrain.mGreenLagrangeStrain[5] ;

	return 0.5*(
			greenLagrangeStrain.mGreenLagrangeStrain[0]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[0]
		   +greenLagrangeStrain.mGreenLagrangeStrain[1]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[1]
		   +greenLagrangeStrain.mGreenLagrangeStrain[2]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[2]
		   +greenLagrangeStrain.mGreenLagrangeStrain[3]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[3]
		   +greenLagrangeStrain.mGreenLagrangeStrain[4]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[4]
		   +greenLagrangeStrain.mGreenLagrangeStrain[5]*secondPiolaKirchhoffStress.mSecondPiolaKirchhoffStress[5]);
}

// calculate coefficients of the material matrix
void NuTo::LinearElastic::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2.*(1.0 + this->mNu));
}

///////////////////////////////////////////////////////////////////////////

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get density
//! @return ... density
double NuTo::LinearElastic::GetDensity() const
{
	return this->mRho;
}

//! @brief ... set density
//! @param rRho ... density
void NuTo::LinearElastic::SetDensity(double rRho)
{
    this->CheckDensity(rRho);
    this->mRho = rRho;
    this->SetParametersValid();
}

//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::LinearElastic::GetYoungsModulus() const
{
	return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::LinearElastic::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::LinearElastic::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::LinearElastic::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::LinearElastic::GetType() const
{
    return NuTo::Constitutive::LINEAR_ELASTIC;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::LinearElastic::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
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

//! @brief ... check if density is positive
//! @param rRho ... density
void NuTo::LinearElastic::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::LinearElastic::CheckDensity] The density must be a positive value.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::LinearElastic::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::LinearElastic::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::LinearElastic::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::LinearElastic::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::LinearElastic::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::LinearElastic::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus: " << this->mE << "\n";
    rLogger << "    Poisson's ratio: " << this->mNu << "\n";
    rLogger << "    Density        : " << this->mRho << "\n";
}

// check parameters
void NuTo::LinearElastic::CheckParameters()const
{
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckDensity(this->mRho);
}

