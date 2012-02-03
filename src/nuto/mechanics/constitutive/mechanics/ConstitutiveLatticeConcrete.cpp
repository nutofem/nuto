// $Id: ConstitutiveLatticeConcrete1D.cpp 112 2009-11-17 16:43:15Z unger3 $

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
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal2x2.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveLatticeConcrete.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain3D.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

NuTo::ConstitutiveLatticeConcrete::ConstitutiveLatticeConcrete() : ConstitutiveLatticeStressStrain()
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
void NuTo::ConstitutiveLatticeConcrete::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
   std::cout << "start serialize ConstitutiveLatticeConcrete" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveLatticeStressStrain)
      & BOOST_SERIALIZATION_NVP(mE)
      & BOOST_SERIALIZATION_NVP(mNu)
      & BOOST_SERIALIZATION_NVP(mRho);
#ifdef DEBUG_SERIALIZATION
   std::cout << "finish serialize ConstitutiveLatticeConcrete" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveLatticeConcrete)
#endif // ENABLE_SERIALIZATION

//  Lattice strain /////////////////////////////////////
//! @brief ... calculate Lattice plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rLatticeStrain ... Lattice strain
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetLatticePlasticStrain(const ElementBase* rElement, int rIp,
								  const LatticeStrain2D& rLatticeStrain, LatticeStrain3D& rLatticePlasticStrain) const
{
	rLatticePlasticStrain.mLatticeStrain[0] = 0.;
	rLatticePlasticStrain.mLatticeStrain[1] = 0.;

	return Error::SUCCESSFUL;
}

//  Lattice strain /////////////////////////////////////
//! @brief ... calculate Lattice plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rLatticeStrain ... Lattice strain
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetLatticePlasticStrain(const ElementBase* rElement, int rIp,
								  const LatticeStrain3D& rLatticeStrain, LatticeStrain3D& rLatticePlasticStrain) const
{
	rLatticePlasticStrain.mLatticeStrain[0] = 0.;
	rLatticePlasticStrain.mLatticeStrain[1] = 0.;
	rLatticePlasticStrain.mLatticeStrain[2] = 0.;

    return Error::SUCCESSFUL;
}

// Lattice stress - Lattice strain /////////////////////////////////////
//! @brief ... calculate Lattice stress from Lattice strain
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
			  const LatticeStrain2D& rLatticeStrain, LatticeStress2D& rLatticeStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetLatticeStressFromLatticeStrain] Check the material parameters.");
    }

    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate Lattice stress
	rLatticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	rLatticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];

    return Error::SUCCESSFUL;
}

// Lattice stress - Lattice strain /////////////////////////////////////
//! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
			  const LatticeStrain2D& rLatticeStrain, LatticeStress3D& rLatticeStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetLatticeStressFromLatticeStrain] Check the material parameters.");
    }

    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate Lattice stress
	rLatticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	rLatticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];
	rLatticeStress.mLatticeStress[2] = 0.;

    return Error::SUCCESSFUL;
}

// Lattice stress - Lattice strain /////////////////////////////////////
//! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
			  const LatticeStrain3D& rLatticeStrain, LatticeStress3D& rLatticeStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetLatticeStressFromLatticeStrain] Check the material parameters.");
    }

    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate Lattice stress
	rLatticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	rLatticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];
	rLatticeStress.mLatticeStress[2] = Et * rLatticeStrain.mLatticeStrain[2];

    return Error::SUCCESSFUL;
}

 //  Damage /////////////////////////////////////
 //! @brief ... calculate isotropic damage from deformation gradient in 2D
 //! @param rElement ... element
 //! @param rIp ... integration point
 //! @param rLatticeStrain ... deformation gradient
 //! @param rDamage ... damage variable
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetDamage(const ElementBase* rElement, int rIp,
                                   const LatticeStrain2D& rLatticeStrain, double& rDamage) const
{
    rDamage=0.;
    return Error::SUCCESSFUL;
}

 //  Damage /////////////////////////////////////
 //! @brief ... calculate isotropic damage from deformation gradient in 3D
 //! @param rElement ... element
 //! @param rIp ... integration point
 //! @param rLatticeStrain ... deformation gradient
 //! @param rDamage ... damage variable
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetDamage(const ElementBase* rElement, int rIp,
                                   const LatticeStrain3D& rLatticeStrain, double& rDamage) const
{
    rDamage=0.;
    return Error::SUCCESSFUL;
}

//! @brief ... calculate the tangent (derivative of the Lattice stresses with respect to the Lattice strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rTangent ... tangent
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetTangent_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain2D& rLatticeStrain,
		ConstitutiveTangentBase* rTangent) const
{
	ConstitutiveTangentLocal2x2 *tangent(rTangent->AsConstitutiveTangentLocal2x2());
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetTangent_LatticeStress_LatticeStrain] Check the material parameters.");
    }

    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // store tangent at the output object
    tangent->mTangent[ 0] = En;
    tangent->mTangent[ 1] = 0;
    tangent->mTangent[ 2] = 0;
    tangent->mTangent[ 3] = Et;

    rTangent->SetSymmetry(true);
    return Error::SUCCESSFUL;
}


//! @brief ... calculate the tangent (derivative of the Lattice stresses with respect to the Lattice strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
//! @param rTangent ... tangent
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetTangent_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain3D& rLatticeStrain,
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
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetTangent_LatticeStress_LatticeStrain] Check the material parameters.");
    }

    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // store tangent at the output object
    tangent->mTangent[ 0] = En;
    tangent->mTangent[ 1] = 0.;
    tangent->mTangent[ 2] = 0.;

    tangent->mTangent[ 3] = 0.;
    tangent->mTangent[ 4] = Et;
    tangent->mTangent[ 5] = 0.;

    tangent->mTangent[ 6] = 0.;
    tangent->mTangent[ 7] = 0.;
    tangent->mTangent[ 8] = Et;

    rTangent->SetSymmetry(true);
    return Error::SUCCESSFUL;
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::UpdateStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
		const LatticeStrain2D& rLatticeStrain) const
{
    //no static data required -> empty routine
    return Error::SUCCESSFUL;
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::UpdateStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
		const LatticeStrain3D& rLatticeStrain) const
{
	//no static data required -> empty routine
    return Error::SUCCESSFUL;
}


//! @brief ... update mp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::UpdateTmpStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
		const LatticeStrain2D& rLatticeStrain) const
{
	//no need to update tmp static data
    return Error::SUCCESSFUL;
}

//! @brief ... update mp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::UpdateTmpStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
		const LatticeStrain3D& rLatticeStrain) const
{
	//no need to update tmp static data
    return Error::SUCCESSFUL;
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveLatticeConcrete::AllocateStaticDataLatticeStress_LatticeStrain2D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
	return 0;
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveLatticeConcrete::AllocateStaticDataLatticeStress_LatticeStrain3D(
		const ElementBase* rElement) const
{
	//no static data return Null-Pointer
	return 0;
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetTotalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain2D& rLatticeStrain, double& rEnergy) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetTotalEnergy_LatticeStress_LatticeStrain] Check the material parameters.");
    }
    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate Lattice stress
    LatticeStress2D latticeStress;
	latticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	latticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];

    rEnergy = 0.5*(
    		rLatticeStrain.mLatticeStrain[0]*latticeStress.mLatticeStress[0]
           +rLatticeStrain.mLatticeStrain[1]*latticeStress.mLatticeStress[1]);

    return Error::SUCCESSFUL;
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetTotalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain3D& rLatticeStrain, double& rEnergy) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetTotalEnergy_LatticeStress_LatticeStrain] Check the material parameters.");
    }
    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate Lattice stress
    LatticeStress2D latticeStress;
	latticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	latticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];
	latticeStress.mLatticeStress[2] = Et * rLatticeStrain.mLatticeStrain[2];

    rEnergy = 0.5*(
    		rLatticeStrain.mLatticeStrain[0]*latticeStress.mLatticeStress[0]
    	   +rLatticeStrain.mLatticeStrain[1]*latticeStress.mLatticeStress[1]
           +rLatticeStrain.mLatticeStrain[2]*latticeStress.mLatticeStress[2]);

    return Error::SUCCESSFUL;
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain2D& rLatticeStrain, double& rEnergy) const
{
	return GetTotalEnergy_LatticeStress_LatticeStrain(rElement, rIp, rLatticeStrain, rEnergy);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeConcrete::GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain3D& rLatticeStrain, double& rEnergy) const
{
	return GetTotalEnergy_LatticeStress_LatticeStrain(rElement, rIp, rLatticeStrain, rEnergy);
}

///////////////////////////////////////////////////////////////////////////

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get density
//! @return ... density
double NuTo::ConstitutiveLatticeConcrete::GetDensity() const
{
	return this->mRho;
}

//! @brief ... set density
//! @param rRho ... density
void NuTo::ConstitutiveLatticeConcrete::SetDensity(double rRho)
{
    this->CheckDensity(rRho);
    this->mRho = rRho;
    this->SetParametersValid();
}

//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::ConstitutiveLatticeConcrete::GetYoungsModulus() const
{
	return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::ConstitutiveLatticeConcrete::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::ConstitutiveLatticeConcrete::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::ConstitutiveLatticeConcrete::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::ConstitutiveLatticeConcrete::GetType() const
{
    return NuTo::Constitutive::LATTICE_CONCRETE;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::ConstitutiveLatticeConcrete::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::LATTICE2D:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is positive
//! @param rRho ... density
void NuTo::ConstitutiveLatticeConcrete::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveLatticeConcrete::CheckDensity] The density must be a positive value.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::ConstitutiveLatticeConcrete::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveLatticeConcrete::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::ConstitutiveLatticeConcrete::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveLatticeConcrete::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveLatticeConcrete::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstitutiveLatticeConcrete::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus: " << this->mE << "\n";
    rLogger << "    Poisson's ratio: " << this->mNu << "\n";
    rLogger << "    Density        : " << this->mRho << "\n";
}

// check parameters
void NuTo::ConstitutiveLatticeConcrete::CheckParameters()const
{
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckDensity(this->mRho);
}

//! @brief ... calculate the elastic parameters for the plane
void NuTo::ConstitutiveLatticeConcrete::CalculateElasticParameters(double& rEn,double& rEt)const
{
    rEn = mE/(1.-2*mNu);
    double alpha = (1.-4.*mNu)/(1.+mNu);
    rEt = rEn*alpha;
}

