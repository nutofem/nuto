

#ifndef CONSTITUTIVESTATICDATAGRADIENTDAMAGE1DFATIGUE_H
#define CONSTITUTIVESTATICDATAGRADIENTDAMAGE1DFATIGUE_H

#include "mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain1D.h"
#include "mechanics/constitutive/mechanics/NonlocalEqStrain.h"


//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Vitaliy Kindrachuk, BAM
//! @date November 2015
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataGradientDamage1DFatigue : public ConstitutiveStaticDataPrevEngineeringStressStrain1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientDamageEngineeringStressFatigue;
public:
	//! @brief constructor
    ConstitutiveStaticDataGradientDamage1DFatigue();

    //! @brief copy constructor
    ConstitutiveStaticDataGradientDamage1DFatigue(ConstitutiveStaticDataGradientDamage1DFatigue const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataGradientDamage1DFatigue* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataGradientDamage1DFatigue& operator= (ConstitutiveStaticDataGradientDamage1DFatigue const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataGradientDamage1DFatigue* AsGradientDamage1DFatigue();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataGradientDamage1DFatigue* AsGradientDamage1DFatigue()const;

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

    //!@brief get prev nonlocal eq strain
    NonlocalEqStrain GetNonlocalEqStrain()
    {
    	return this->mPrevNonlocalEqStrain;
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

    //! @brief previous nonlocal eq strain
    NuTo::NonlocalEqStrain mPrevNonlocalEqStrain;

};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataGradientDamage1DFatigue)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAGRADIENTDAMAGE1DFATIGUE_H
