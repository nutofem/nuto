#pragma once

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Thomas Titscher, BAM
//! @date March 2016
namespace NuTo
{
class MisesPlasticityEngineeringStress;
class IpDataStaticDataBase;

template <int TDim>
class ConstitutiveStaticDataMisesPlasticity : virtual public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
    friend class MisesPlasticityEngineeringStress;
public:
	//! @brief constructor
	ConstitutiveStaticDataMisesPlasticity();

    //! @brief copy constructor
	ConstitutiveStaticDataMisesPlasticity(ConstitutiveStaticDataMisesPlasticity const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    virtual ConstitutiveStaticDataMisesPlasticity* Clone()const
    {
    	return new ConstitutiveStaticDataMisesPlasticity(*this);
    }

    //! @brief check, if the static data is compatible with a givenr element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //! @brief assignment operator
    ConstitutiveStaticDataMisesPlasticity& operator= (ConstitutiveStaticDataMisesPlasticity const& rOther) = default;


          NuTo::ConstitutiveStaticDataMisesPlasticity<1>* AsConstitutiveStaticDataMisesPlasticity1D() override
        { throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong dimension"); }
    const NuTo::ConstitutiveStaticDataMisesPlasticity<1>* AsConstitutiveStaticDataMisesPlasticity1D() const override
        { throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong dimension"); }

          NuTo::ConstitutiveStaticDataMisesPlasticity<2>* AsConstitutiveStaticDataMisesPlasticity2D() override
        { throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong dimension"); }
    const NuTo::ConstitutiveStaticDataMisesPlasticity<2>* AsConstitutiveStaticDataMisesPlasticity2D() const override
        { throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong dimension"); }

         NuTo::ConstitutiveStaticDataMisesPlasticity<3>* AsConstitutiveStaticDataMisesPlasticity3D() override
        { throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong dimension"); }
    const NuTo::ConstitutiveStaticDataMisesPlasticity<3>* AsConstitutiveStaticDataMisesPlasticity3D() const override
        { throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong dimension"); }

protected:
    //! @brief accumulated plastic strain (is not always equivalent to epsilon_p)
    double mEpsilonPEq;

    //! @brief plastic strain
    EngineeringStrain<TDim> mEpsilonP;

    //! @brief back stress
    EngineeringStress<TDim> mSigmaB;
};

}

