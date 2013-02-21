// $Id: ConstitutiveStaticDataGradientDamagePlasticity1D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAGRADIENTDAMAGEPLASTICITY1D_H
#define CONSTITUTIVESTATICDATAGRADIENTDAMAGEPLASTICITY1D_H

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain1D.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataGradientDamagePlasticity1D : public ConstitutiveStaticDataPrevEngineeringStressStrain1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientDamagePlasticityEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataGradientDamagePlasticity1D();

    //! @brief copy constructor
    ConstitutiveStaticDataGradientDamagePlasticity1D(ConstitutiveStaticDataGradientDamagePlasticity1D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataGradientDamagePlasticity1D* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataGradientDamagePlasticity1D& operator= (ConstitutiveStaticDataGradientDamagePlasticity1D const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataGradientDamagePlasticity1D* AsGradientDamagePlasticity1D();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataGradientDamagePlasticity1D* AsGradientDamagePlasticity1D()const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief local accumulated plastic strain
    double mKappa;

    //! @brief previous (after update) total strain component in radial direction (e22=e33)
    double mEpsilonTotRadial;

    //! @brief plastic strain (e11 and the radial component of the plastic strain e22=e33)
    double mEpsilonP[2];
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAGRADIENTDAMAGEPLASTICITY1D_H
