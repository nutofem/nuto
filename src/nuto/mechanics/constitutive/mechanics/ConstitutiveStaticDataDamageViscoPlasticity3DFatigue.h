// $Id: ConstitutiveStaticDataDamageViscoPlasticity3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATADAMAGEVISCOPLASTICITY3DFATIGUE_H
#define CONSTITUTIVESTATICDATADAMAGEVISCOPLASTICITY3DFATIGUE_H

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3D.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class DamageViscoPlasticityHardeningEngineeringStress;
class DamageViscoPlasticityEngineeringStress;
class IpDataStaticDataBase;

class ConstitutiveStaticDataDamageViscoPlasticity3DFatigue : public ConstitutiveStaticDataDamageViscoPlasticity3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class DamageViscoPlasticityEngineeringStress;
    friend class DamageViscoPlasticityHardeningEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataDamageViscoPlasticity3DFatigue();

    //! @brief copy constructor
    ConstitutiveStaticDataDamageViscoPlasticity3DFatigue(ConstitutiveStaticDataDamageViscoPlasticity3DFatigue const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* Clone()const
    {
    	return new ConstitutiveStaticDataDamageViscoPlasticity3DFatigue(*this);
    }

    //! @brief assignment operator
    ConstitutiveStaticDataDamageViscoPlasticity3DFatigue& operator= (ConstitutiveStaticDataDamageViscoPlasticity3DFatigue const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as damage viscoplasticity static data
    ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* AsDamageViscoPlasticity3DFatigue();

    //!@ brief reinterpret as damage viscoplasticity static data
    const ConstitutiveStaticDataDamageViscoPlasticity3DFatigue* AsDamageViscoPlasticity3DFatigue()const;

    //!@brief get mKappaInelasticFatigue
    double GetKappaInelasticFatigue()
    {
    	return this->mKappaInelasticFatigue;
    }

    //!@brief set mKappaInelasticFatigue
    void SetKappaInelasticFatigue(double KappaInelastic)
    {
    	this->mKappaInelasticFatigue = KappaInelastic;
    }

    //!@brief get mOmegaComprFatigue
    double GetOmegaComprFatigue()
    {
    	return this->mOmegaComprFatigue;
    }

    //!@brief set mOmegaComprFatigue
    void SetOmegaComprFatigue(double OmegaCompr)
    {
    	this->mOmegaComprFatigue = OmegaCompr;
    }

    //!@brief get mPrevHardeningFatigue
    double GetPrevHardeningFatigue()
    {
    	return this->mPrevHardeningFatigue;
    }

    //!@brief set mPrevHardeningFatigue
    void SetPrevHardeningFatigue(double PrevHardening)
    {
    	this->mPrevHardeningFatigue = PrevHardening;
    }

    //!@brief get mEpsilonVpFatigue
    EngineeringStrain3D GetEpsilonVpFatigue()
    {
    	return this->mEpsilonVpFatigue;
    }

    //!@brief set mEpsilonVpFatigue
    void SetEpsilonVpFatigue(EngineeringStrain3D EpsilonVp)
    {
    	this->mEpsilonVpFatigue = EpsilonVp;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief distributed accumulated inelastic strain
    double mKappaInelasticFatigue;

    //! @brief local damage variable associated with plastic strain (compressive damage)
    double mOmegaComprFatigue;

    //! @brief plastic strain and previous strain
    EngineeringStrain3D mEpsilonPFatigue;

    //! @brief viscoplastic strain
    EngineeringStrain3D mEpsilonVpFatigue;

    //! @brief previous strain
    EngineeringStrain3D mPrevStrainFatigue;

    //! @brief previous stress
    EngineeringStress3D mPrevStressFatigue;

    //! @brief hardening state variable
    double mPrevHardeningFatigue;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataDamageViscoPlasticity3DFatigue)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATADAMAGEVISCOPLASTICITY3DFATIGUE_H
