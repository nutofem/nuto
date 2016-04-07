#pragma once

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Thomas Titscher, BAM
//! @date November 2014
namespace NuTo
{

class ConstitutiveStaticDataGradientDamage : public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientDamageEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataGradientDamage();

    //! @brief copy constructor
    ConstitutiveStaticDataGradientDamage(ConstitutiveStaticDataGradientDamage const& rOther) = default;

    //! @brief clones (copies) the data
    ConstitutiveStaticDataGradientDamage* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataGradientDamage& operator= (ConstitutiveStaticDataGradientDamage const& rOther) = default;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataGradientDamage* AsGradientDamage();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataGradientDamage* AsGradientDamage()const;

    void SetKappa(double rKappa)
    {
        mKappa = rKappa;
    }

    double GetKappa() const
    {
        return mKappa;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION



protected:
    //! @brief maximal value of nonlocal eq strain
    double mKappa;

};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataGradientDamage)
#endif // ENABLE_SERIALIZATION
