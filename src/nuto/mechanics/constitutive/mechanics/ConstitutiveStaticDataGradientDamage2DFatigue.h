

#ifndef CONSTITUTIVESTATICDATAGRADIENTDAMAGE2DFATIGUE_H
#define CONSTITUTIVESTATICDATAGRADIENTDAMAGE2DFATIGUE_H

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"


//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Vitaliy Kindrachuk, BAM
//! @date November 2015
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataGradientDamage2DFatigue : public ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientDamageEngineeringStressFatigue;
public:
	//! @brief constructor
    ConstitutiveStaticDataGradientDamage2DFatigue();

    //! @brief copy constructor
    ConstitutiveStaticDataGradientDamage2DFatigue(ConstitutiveStaticDataGradientDamage2DFatigue const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataGradientDamage2DFatigue* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataGradientDamage2DFatigue& operator= (ConstitutiveStaticDataGradientDamage2DFatigue const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataGradientDamage2DFatigue* AsGradientDamage2DFatigue();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataGradientDamage2DFatigue* AsGradientDamage2DFatigue()const;

    //!@brief get mKappaFatigue
    double GetKappaFatigue()
    {
    	return this->mKappaFatigue;
    }

    //!@brief get mOmegaFatigue
    double GetOmegaFatigue()
    {
    	return this->mOmegaFatigue;
    }

    //!@brief get mKappa
    double GetKappa()
    {
    	return this->mKappa;
    }

    //!@brief get mOmega
    double GetOmega()
    {
    	return this->mOmega;
    }

    //!@brief get mPrevNonlocalEqStrain
    double GetPrevNonlocalEqStrain()
    {
    	return this->mPrevNonlocalEqStrain.GetValue(0);
    }

    //!@brief get mPrevNonlocalEqStrainFatigue
    double GetPrevNonlocalEqStrainFatigue()
    {
    	return this->mPrevNonlocalEqStrainFatigue.GetValue(0);
    }

    //!@brief get mPrevStrainFatigue
    EngineeringStrain2D GetPrevStrainFatigue()
    {
    	return this->mPrevStrainFatigue;
    }

    //!@brief get mPrevStressFatigue
    EngineeringStress2D GetPrevStressFatigue()
    {
    	return this->mPrevStressFatigue;
    }

    void SetKappaFatigue(double rKappa)
    {
        mKappaFatigue = rKappa;
    }

    void SetOmegaFatigue(double rOmega)
    {
        mOmegaFatigue = rOmega;
    }

    void SetKappa(double rKappa)
    {
        mKappa = rKappa;
    }

    void SetOmega(double rOmega)
    {
        mOmega = rOmega;
    }

    //!@brief save static data to their relevant fatigue counterparts
    void FatigueSaveStaticData()
    {
    	this->SetOmegaFatigue(this->mOmega);
    	this->SetKappaFatigue(this->mKappa);
    	this->mPrevStrainFatigue = this->mPrevStrain;
    	this->mPrevStressFatigue = this->mPrevSigma;
    	this->mPrevNonlocalEqStrainFatigue = this->mPrevNonlocalEqStrain;
    }

    //!@brief save static data to their relevant fatigue counterparts
    void FatigueRestoreStaticData()
    {
    	this->SetOmega(this->mOmegaFatigue);
    	this->SetKappa(this->mKappaFatigue);
    	this->mPrevStrain = this->mPrevStrainFatigue;
    	this->mPrevSigma = this->mPrevStressFatigue;
    	this->mPrevNonlocalEqStrain = this->mPrevNonlocalEqStrainFatigue;
    }

    //!@brief extrapolate static data except of mPrevSigma and mPrevStrain
    //!@brief mPrevSigma and mPrevStrain should be calculated after finding the equilibrium with the extrapolated static data
    void FatigueExtrapolateStaticData(int rNumber)
	{
    	// the kappa should be extrapolated, damage should be calculated from kappa in the IP->Evaluate
    	double DeltaKappa;

    	// calculate cyclic change of kappa
    	DeltaKappa = this->mKappa - this->mKappaFatigue;

    	// linear extrapolation of kappa
    	this->mKappa += rNumber*DeltaKappa;
	}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION



protected:
    //! @brief damage variable
    double mOmega;

    //! @brief maximal value of nonlocal eq strain
    double mKappa;

    //! @brief fatigue damage variable
    double mOmegaFatigue;


    //! @brief maximal value of nonlocal eq strain
    double mKappaFatigue;

    //! @brief previous strain
    EngineeringStrain2D mPrevStrainFatigue;

    //! @brief previous stress
    EngineeringStress2D mPrevStressFatigue;

    //! @brief previous nonlocal eq strain
    NuTo::NonlocalEqStrain mPrevNonlocalEqStrain;

    //! @brief previous nonlocal eq strain
    NuTo::NonlocalEqStrain mPrevNonlocalEqStrainFatigue;

};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataGradientDamage2DFatigue)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAGRADIENTDAMAGE2DFATIGUE_H
