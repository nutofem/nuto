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
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataLatticeConcrete2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain3D.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

NuTo::ConstitutiveLatticeConcrete::ConstitutiveLatticeConcrete() : ConstitutiveLatticeStressStrain()
{
	mE = 30000.;
	mNu = 0.;
	mRho = 0.;
	mAlpha = 0.25;
	mSigmaT = 3.0;
	mSigmaN = 3.0;
	mG = 0.1;
	mNt = 0.2;
	mMu = 0.0;
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
      & BOOST_SERIALIZATION_NVP(mRho)
      & BOOST_SERIALIZATION_NVP(mAlpha)
      & BOOST_SERIALIZATION_NVP(mSigmaT)
      & BOOST_SERIALIZATION_NVP(mSigmaN)
      & BOOST_SERIALIZATION_NVP(mG);
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
#define tol 1e-14
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

	//calculate the boundary stress
	if (rLatticeStrain.mLatticeStrain[0]>0)
	{
	    //tension cutoff
		CutOffTension2D(rElement, rIp, rLatticeStrain, &rLatticeStress, 0);
	}
	else
	{
	    //compression cutoff
		CutOffCompression2D(rElement, rIp, rLatticeStrain, &rLatticeStress, 0);
	}

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
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::GetTangent_LatticeStress_LatticeStrain] Check the material parameters.");
    }

	//calculate the boundary stress
	if (rLatticeStrain.mLatticeStrain[0]>0)
	{
	    //tension cutoff
		CutOffTension2D(rElement, rIp, rLatticeStrain, 0, rTangent->AsConstitutiveTangentLocal2x2());
	}
	else
	{
	    //compression cutoff
		CutOffCompression2D(rElement, rIp, rLatticeStrain, 0, rTangent->AsConstitutiveTangentLocal2x2());
	}

    rTangent->SetSymmetry(false);
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
	double epsilonE;
	if (rLatticeStrain.mLatticeStrain[0]>0)
	{
		//tension
		epsilonE = sqrt(rLatticeStrain.mLatticeStrain[0]*rLatticeStrain.mLatticeStrain[0]+mAlpha*rLatticeStrain.mLatticeStrain[1]*rLatticeStrain.mLatticeStrain[1]);
	}
	else
	{
		//compression
		epsilonE = sqrt(mAlpha)*fabs(rLatticeStrain.mLatticeStrain[1]);
	}

	ConstitutiveStaticDataLatticeConcrete2D*
	    staticDataPtr(rElement->GetStaticData(rIp)->AsConstitutiveStaticDataLatticeConcrete2D());
	if (epsilonE>staticDataPtr->mEpsilonMax)
	{
		staticDataPtr->mEpsilonMax = epsilonE;
	}

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
	return new NuTo::ConstitutiveStaticDataLatticeConcrete2D();
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveLatticeConcrete::AllocateStaticDataLatticeStress_LatticeStrain3D(
		const ElementBase* rElement) const
{
	return 0;
	//return new NuTo::ConstitutiveStaticDataLatticeConcrete3D();
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
/*    rEn = mE/(1.-2*mNu);
    double alpha = (1.-4.*mNu)/(1.+mNu);
    rEt = rEn*alpha;
*/
	rEn = mE*(4.+mAlpha)/(2.+3.*mAlpha);
	rEt = rEn*mAlpha;
}

/*//! @brief ... calculate the boundary stress in tension
void NuTo::ConstitutiveLatticeConcrete::CutOffTension2D(const ElementBase* rElement, int rIp,
		  const LatticeStrain2D& rLatticeStrain, LatticeStress2D&rLatticeStress)const
{
    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate trial lattice stress
	rLatticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	rLatticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];

	//get the static data
	const ConstitutiveStaticDataLatticeConcrete2D *oldStaticData = (rElement->GetStaticData(rIp))->AsConstitutiveStaticDataLatticeConcrete2D();
	double epsilonMax(oldStaticData->mEpsilonMax);
	double ledge(rElement->GetIpEdgeLength(rIp));

	//calculate omega
	double omega;
	if (fabs(rLatticeStrain.mLatticeStrain[1])>tol)
	{
		omega = atan(rLatticeStrain.mLatticeStrain[0]/(sqrt(mAlpha)*fabs(rLatticeStrain.mLatticeStrain[1])));
	}
	else
	{
		omega = 0;
	}

	//calculate maximum effective strain
	double epsilonE = sqrt(rLatticeStrain.mLatticeStrain[0]*rLatticeStrain.mLatticeStrain[0]+mAlpha*rLatticeStrain.mLatticeStrain[1]*rLatticeStrain.mLatticeStrain[1]);
	if (epsilonE>epsilonMax-toleranceEpsilonMax)
		epsilonMax=epsilonE;

	//calculate sigma0
	double rst = mSigmaT/mSigmaN;
	double sinOmega = sin(omega);
	double cosOmega = cos(omega);

	double sigma0(mSigmaN*(sqrt(rst*rst*cosOmega*cosOmega/mAlpha+sinOmega*sinOmega)));

	double epsilon0 = sigma0/En;


	double EmaxMinusE0 = epsilonMax-epsilon0;

	double sigmaBT;
	if (EmaxMinusE0>0)
	{
		double lt = 2.*En*mG/(mSigmaN*mSigmaN);
		double Ht = 2.*En/(lt/ledge-1.);

		double H0 = Ht*pow(2.*omega/M_PI,mNt);
		sigmaBT = sigma0*exp(-H0*(EmaxMinusE0)/sigma0);
	}
	else
	{
		sigmaBT = sigma0;
	}

	double sigmaE = sqrt(rLatticeStress.mLatticeStress[0]*rLatticeStress.mLatticeStress[0]+rLatticeStress.mLatticeStress[1]*rLatticeStress.mLatticeStress[1]/mAlpha);
	if (sigmaE>sigmaBT)
	{
		rLatticeStress.mLatticeStress[0] = sigmaBT * rLatticeStrain.mLatticeStrain[0]/epsilonE;
		rLatticeStress.mLatticeStress[1] = sigmaBT * mAlpha * rLatticeStrain.mLatticeStrain[1]/epsilonE;
	}

}
*/


//#define SHOWSTRESS
//! @brief ... calculate the stiffness in tension
void NuTo::ConstitutiveLatticeConcrete::CutOffTension2D(const ElementBase* rElement, int rIp,
		  const LatticeStrain2D& rLatticeStrain, LatticeStress2D *rLatticeStress, ConstitutiveTangentLocal2x2 *rTangent)const
{
    assert(rLatticeStrain.mLatticeStrain[0]>=0);

	//Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);
    //LatticeStress2D latticeStress;


	LatticeStrain2D myLatticeStrain(rLatticeStrain);

	//get the static data
	const ConstitutiveStaticDataLatticeConcrete2D *oldStaticData = (rElement->GetStaticData(rIp))->AsConstitutiveStaticDataLatticeConcrete2D();
	double ledge(rElement->GetIpEdgeLength(rIp));

	//calculate maximum effective strain
	double epsilonE = sqrt(rLatticeStrain.mLatticeStrain[0]*rLatticeStrain.mLatticeStrain[0]+mAlpha*rLatticeStrain.mLatticeStrain[1]*rLatticeStrain.mLatticeStrain[1]);

	bool isLoading;
	double epsilonMax(oldStaticData->mEpsilonMax);
	if (epsilonE>epsilonMax-tol)
	{
		epsilonMax=epsilonE;
		isLoading = true;
	}
	else
		isLoading = false;

	//calculate omega
	double omega;
	double dOmegadEpsilon[2];
	double dEpsilonEdEpsilon[2];
	double signTangentStrain(rLatticeStrain.mLatticeStrain[1]<0 ? -1 : 1);
	if (fabs(rLatticeStrain.mLatticeStrain[1])>tol)
	{
		//tension with shear
		omega = atan(rLatticeStrain.mLatticeStrain[0]/(sqrt(mAlpha)*fabs(rLatticeStrain.mLatticeStrain[1])));
		dOmegadEpsilon[0] =  signTangentStrain*rLatticeStrain.mLatticeStrain[1]*sqrt(mAlpha)/(epsilonE* epsilonE);
		dOmegadEpsilon[1] = -rLatticeStrain.mLatticeStrain[0]*sqrt(mAlpha)/(epsilonE* epsilonE)*signTangentStrain;

		dEpsilonEdEpsilon[0] = rLatticeStrain.mLatticeStrain[0]/epsilonE;
		dEpsilonEdEpsilon[1] = mAlpha*rLatticeStrain.mLatticeStrain[1]/epsilonE;
	}
	else
	{
		//pure tension
		omega = M_PI*0.5;

		dOmegadEpsilon[0] = 0.;
		dOmegadEpsilon[1] = 0.;

		if (epsilonE<tol)
		{
			//this is eventually not defined, but for epsilonE==0 we should not reach the boundary
			dEpsilonEdEpsilon[0] = 0.;
			dEpsilonEdEpsilon[1] = 0;
		}
		else
		{
			dEpsilonEdEpsilon[0] = rLatticeStrain.mLatticeStrain[0]/epsilonE;
			dEpsilonEdEpsilon[1] = mAlpha*rLatticeStrain.mLatticeStrain[1]/epsilonE;
		}
	}

	//check dOmegadEpsilon
//	double dOmegadEpsilonCDF[2];
//	double delta(1e-9);
//	for (int count=0; count<2; count++)
//	{
//		myLatticeStrain.mLatticeStrain[count]+=delta;
//		double omega2 = atan(myLatticeStrain.mLatticeStrain[0]/(sqrt(mAlpha)*fabs(myLatticeStrain.mLatticeStrain[1])));
//		myLatticeStrain.mLatticeStrain[count]-=delta;
//		dOmegadEpsilonCDF[count] = (omega2-omega)/delta;
//	}
//	std::cout << "dOmegadEpsilon ana " << dOmegadEpsilon[0] << " " << dOmegadEpsilon[1] << "\n";
//	std::cout << "dOmegadEpsilon cdf " << dOmegadEpsilonCDF[0] << " " << dOmegadEpsilonCDF[1] << "\n";


	//calculate sigma0
	double rst = mSigmaT/mSigmaN;
	double sinOmega = sin(omega);
	double cosOmega = cos(omega);

	double sigma0(mSigmaN*(sqrt(rst*rst*cosOmega*cosOmega/mAlpha+sinOmega*sinOmega)));
	double dSigma0dOmega (mSigmaN*mSigmaN/sigma0*cosOmega*sinOmega*(1-rst*rst/mAlpha));

//	double sigma0b (mSigmaN*(sqrt(mAlpha*rst*rst*cos(omega+delta)*cos(omega+delta)+sin(omega+delta)*sin(omega+delta))));
//	std::cout << "dSigma0dOmega ana: " << dSigma0dOmega << " cdf: " << (sigma0b-sigma0)/delta << "\n";

	double epsilon0 = sigma0/En;

	if (epsilonMax<epsilon0)
	{

	    //elastic loading
	    if (rLatticeStress!=0)
	    {
	    	rLatticeStress->mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	    	rLatticeStress->mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];
	    }
	    if (rTangent!=0)
	    {
			rTangent->mTangent[0] = En;
			rTangent->mTangent[1] = 0.;
			rTangent->mTangent[2] = 0.;
			rTangent->mTangent[3] = Et;
	    }
#ifdef SHOWSTRESS
	    std::cout << rElement->ElementGetId() << ":" << rIp << " tension elastic" << "\n";
#endif
	}
	else
	{

		double EmaxMinusE0 = epsilonMax-epsilon0;

		double sigmaBT;
		double dSigmaBTdEpsilon[2];
		if (EmaxMinusE0>0)
		{
			double lt = 2.*En*mG/(mSigmaN*mSigmaN);
			double Ht = 2.*En/(lt/ledge-1.);
			if (lt/ledge<1.1)
				throw MechanicsException("[NuTo::ConstitutiveLatticeConcrete::CutOffTension2D] snap back on material level - decrease your element size.");
			double H0 = Ht*pow(2.*omega/M_PI,mNt);
			double expValue = exp(-H0*(EmaxMinusE0)/sigma0);
			sigmaBT = sigma0*expValue;

			double dH0dOmega;
			if (fabs(omega)>tol)
			{
				dH0dOmega = mNt/omega*H0;
	//std::cout << "dH0dOmega ana : " << dH0dOmega << " cdf " <<  (Ht*pow(2.*(omega+delta)/M_PI,mNt)-H0)/delta << "\n";
			}
			else
			{
				dH0dOmega = Ht*pow(2*tol/M_PI,mNt-1.)*2./M_PI;
			}

			double dExpValuedOmega= expValue*(H0*epsilonMax*dSigma0dOmega/(sigma0*sigma0)-dH0dOmega*EmaxMinusE0/sigma0);

	//		{
	//			double sigma02 = mSigmaN*(sqrt(mAlpha*rst*rst*cos(omega+delta)*cos(omega+delta)+sin(omega+delta)*sin(omega+delta)));
	//			double H02 = Ht*pow(2.*(omega+delta)/M_PI,mNt);
	//			double expValue2 = exp(-H02*(epsilonMax/sigma02-1./En));
	//			std::cout << "dExpValuedOmega ana : " << dExpValuedOmega << " cdf : " << (expValue2-expValue)/delta << "\n";
	//		}

			dSigmaBTdEpsilon[0] = (dSigma0dOmega*expValue+sigma0*dExpValuedOmega)*dOmegadEpsilon[0];
			dSigmaBTdEpsilon[1] = (dSigma0dOmega*expValue+sigma0*dExpValuedOmega)*dOmegadEpsilon[1];
			if (isLoading)
			{
				double dExpvaluedEpsilonMax = -expValue*H0/sigma0;
	//			{
	//				double EmaxMinusE02 = (epsilonMax+delta)-epsilon0;
	//				double expValue2 = exp(-H0*(EmaxMinusE02)/sigma0);
	//				std::cout << "dExpvaluedEpsilonMax ana : " << dExpvaluedEpsilonMax << " cdf : " << (expValue2-expValue)/delta << "\n";
	//			}

				dSigmaBTdEpsilon[0] += sigma0*dExpvaluedEpsilonMax*dEpsilonEdEpsilon[0];
				dSigmaBTdEpsilon[1] += sigma0*dExpvaluedEpsilonMax*dEpsilonEdEpsilon[1];
			}
		}
		else
		{
			sigmaBT = sigma0;
			dSigmaBTdEpsilon[0] = dSigma0dOmega*dOmegadEpsilon[0];
			dSigmaBTdEpsilon[1] = dSigma0dOmega*dOmegadEpsilon[1];
		}

	//	double dSigmaBTdEpsilon2[2];
	//	for (int count=0; count<2; count ++)
	//	{
	//		LatticeStrain2D latticeStrain2(rLatticeStrain);
	//		LatticeStress2D latticeStress2;
	//		latticeStrain2.mLatticeStrain[count] +=delta;
	//		double sigmaBT2;
	//
	//		CutOffTension2D(rElement, rIp,latticeStrain2, latticeStress2, sigmaBT2);
	//		dSigmaBTdEpsilon2[count] = (sigmaBT2-sigmaBT)/delta;

	//	}
	//	std::cout << "dSigmaBTdEpsilon ana: " << dSigmaBTdEpsilon[0] << " " << dSigmaBTdEpsilon[1] ;
	//    std::cout << " cdf " << dSigmaBTdEpsilon2[0] << " " << dSigmaBTdEpsilon2[1] << "\n";

	//	double sigmaE = sqrt(latticeStress.mLatticeStress[0]*latticeStress.mLatticeStress[0]+latticeStress.mLatticeStress[1]*latticeStress.mLatticeStress[1]/mAlpha);
		if (isLoading)
		{
			assert(fabs(epsilonE)>tol);

			if (rLatticeStress!=0)
			{
				rLatticeStress->mLatticeStress[0] = sigmaBT * rLatticeStrain.mLatticeStrain[0]/epsilonE;
				rLatticeStress->mLatticeStress[1] = sigmaBT * mAlpha * rLatticeStrain.mLatticeStrain[1]/epsilonE;
			}

			if (rTangent!=0)
			{
				rTangent->mTangent[0] = dSigmaBTdEpsilon[0]*rLatticeStrain.mLatticeStrain[0]/epsilonE+
						sigmaBT*(epsilonE-rLatticeStrain.mLatticeStrain[0]*rLatticeStrain.mLatticeStrain[0]/epsilonE)/(epsilonE*epsilonE);
				rTangent->mTangent[2] = dSigmaBTdEpsilon[1]*rLatticeStrain.mLatticeStrain[0]/epsilonE-
						sigmaBT*(rLatticeStrain.mLatticeStrain[0]*rLatticeStrain.mLatticeStrain[1]*mAlpha)/(epsilonE*epsilonE*epsilonE);
				rTangent->mTangent[1] = dSigmaBTdEpsilon[0]*mAlpha*rLatticeStrain.mLatticeStrain[1]/epsilonE-
						sigmaBT*(mAlpha*rLatticeStrain.mLatticeStrain[0]*rLatticeStrain.mLatticeStrain[1])/(epsilonE*epsilonE*epsilonE);
				rTangent->mTangent[3] = dSigmaBTdEpsilon[1]*rLatticeStrain.mLatticeStrain[1]*mAlpha/epsilonE+
						sigmaBT*mAlpha*(epsilonE-mAlpha*rLatticeStrain.mLatticeStrain[1]*rLatticeStrain.mLatticeStrain[1]/epsilonE)/(epsilonE*epsilonE);
			}
#ifdef SHOWSTRESS
			std::cout << rElement->ElementGetId() << ":" << rIp << " tension loading" << "\n";
#endif
		}
		else
		{
			//unloading back to the origin
			double help = sigmaBT/epsilonMax;

			if (rLatticeStress->mLatticeStress!=0)
			{
				rLatticeStress->mLatticeStress[0] = help * rLatticeStrain.mLatticeStrain[0];
				rLatticeStress->mLatticeStress[1] = help * mAlpha * rLatticeStrain.mLatticeStrain[1];
			}


			if (rTangent!=0)
			{
				rTangent->mTangent[0] = dSigmaBTdEpsilon[0]*rLatticeStrain.mLatticeStrain[0]/epsilonMax + help;
				rTangent->mTangent[1] = dSigmaBTdEpsilon[0]*mAlpha*rLatticeStrain.mLatticeStrain[1]/epsilonMax;
				rTangent->mTangent[2] = dSigmaBTdEpsilon[1]*rLatticeStrain.mLatticeStrain[0]/epsilonMax;
				rTangent->mTangent[3] = dSigmaBTdEpsilon[1]*mAlpha*rLatticeStrain.mLatticeStrain[1]/epsilonMax + help*mAlpha;
			}
#ifdef SHOWSTRESS
			std::cout << rElement->ElementGetId() << ":" << rIp << " tension unloading" << "\n";
#endif
		}
	}
#ifdef SHOWSTRESS
	if (rLatticeStress!=0)
	    std::cout << "stress " << rLatticeStress->mLatticeStress[0] << " " << rLatticeStress->mLatticeStress[1] << "\n";
#endif
}

//! @brief ... calculate the stiffness in tension
void NuTo::ConstitutiveLatticeConcrete::CutOffCompression2D(const ElementBase* rElement, int rIp,
		  const LatticeStrain2D& rLatticeStrain, LatticeStress2D *rLatticeStress, ConstitutiveTangentLocal2x2 *rTangent)const
{
	assert(rLatticeStrain.mLatticeStrain[0]<=0);
    //Calculate parameters
    double En,Et;
    CalculateElasticParameters(En,Et);

    // calculate trial lattice stress
    //LatticeStress2D latticeStress;
    //latticeStress.mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
    //latticeStress.mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];

	//get the static data
	const ConstitutiveStaticDataLatticeConcrete2D *oldStaticData = (rElement->GetStaticData(rIp))->AsConstitutiveStaticDataLatticeConcrete2D();
	bool isLoading;
	double epsilonMax(oldStaticData->mEpsilonMax);
	double epsilonE(fabs(rLatticeStrain.mLatticeStrain[1])*sqrt(mAlpha));
	if (epsilonE>epsilonMax-tol)
	{
		epsilonMax=epsilonE;
		isLoading = true;
	}
	else
		isLoading = false;

	//cohesive plus friction part

	//calculate sqrt Alpha first
	double sigma0  = (mSigmaT - mMu * En * rLatticeStrain.mLatticeStrain[0])/sqrt(mAlpha);
	double epsilon0 = sigma0/En;

	if (epsilonMax<epsilon0)
	{
	    //elastic loading
	    if (rLatticeStress!=0)
	    {
	    	rLatticeStress->mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
	    	rLatticeStress->mLatticeStress[1] = Et * rLatticeStrain.mLatticeStrain[1];
	    }
	    if (rTangent!=0)
	    {
			rTangent->mTangent[0] = En;
			rTangent->mTangent[1] = 0.;
			rTangent->mTangent[2] = 0.;
			rTangent->mTangent[3] = Et;
	    }
#ifdef SHOWSTRESS
	    std::cout << rElement->ElementGetId() << ":" << rIp << " compression elastic" << "\n";
#endif
	}
	else
	{

		if (isLoading)
		{
			if (rLatticeStress!=0)
			{
				rLatticeStress->mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
				rLatticeStress->mLatticeStress[1] = sigma0*sqrt(mAlpha)*(rLatticeStrain.mLatticeStrain[1]<0 ? -1 : 1);
			}

			if (rTangent!=0)
			{
				rTangent->mTangent[0] = En;
				rTangent->mTangent[2] = 0;
				rTangent->mTangent[1] = -mMu*En*(rLatticeStrain.mLatticeStrain[1]<0 ? -1 : 1);
				rTangent->mTangent[3] = 0;
			}
#ifdef SHOWSTRESS
			std::cout << rElement->ElementGetId() << ":" << rIp << " compression loading" << "\n";
#endif
		}
		else
		{
			//unloading back to the origin, only for the tangential part
			if (rLatticeStress->mLatticeStress!=0)
			{
				rLatticeStress->mLatticeStress[0] = En * rLatticeStrain.mLatticeStrain[0];
				rLatticeStress->mLatticeStress[1] = sigma0/epsilonMax * mAlpha * rLatticeStrain.mLatticeStrain[1];
			}


			if (rTangent!=0)
			{
				rTangent->mTangent[0] = En;
				rTangent->mTangent[1] = - mMu * En*sqrt(mAlpha)/epsilonMax * rLatticeStrain.mLatticeStrain[1];
				rTangent->mTangent[2] = 0;
				rTangent->mTangent[3] = sigma0/epsilonMax * mAlpha;
			}
#ifdef SHOWSTRESS
			std::cout << rElement->ElementGetId() << ":" << rIp << " compression unloading" << "\n";
#endif
		}
	}
#ifdef SHOWSTRESS
	if (rLatticeStress!=0)
		std::cout << "stress " << rLatticeStress->mLatticeStress[0] << " " << rLatticeStress->mLatticeStress[1] << "\n";
#endif
}
