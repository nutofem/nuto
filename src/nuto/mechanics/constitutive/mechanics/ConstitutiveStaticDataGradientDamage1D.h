

#ifndef CONSTITUTIVESTATICDATAGRADIENTDAMAGE1D_H
#define CONSTITUTIVESTATICDATAGRADIENTDAMAGE1D_H

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Thomas Titscher, BAM
//! @date November 2014
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataGradientDamage1D : public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientDamageEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataGradientDamage1D();

    //! @brief copy constructor
    ConstitutiveStaticDataGradientDamage1D(ConstitutiveStaticDataGradientDamage1D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataGradientDamage1D* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataGradientDamage1D& operator= (ConstitutiveStaticDataGradientDamage1D const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataGradientDamage1D* AsGradientDamage1D();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataGradientDamage1D* AsGradientDamage1D()const;


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
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataGradientDamage1D)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAGRADIENTDAMAGE1D_H
