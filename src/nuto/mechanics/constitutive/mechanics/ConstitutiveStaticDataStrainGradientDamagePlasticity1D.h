// $Id: ConstitutiveStaticDataGradientDamagePlasticity1D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATASTRAINGRADIENTDAMAGEPLASTICITY1D_H
#define CONSTITUTIVESTATICDATASTRAINGRADIENTDAMAGEPLASTICITY1D_H

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain1D.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataStrainGradientDamagePlasticity1D : public ConstitutiveStaticDataPrevEngineeringStressStrain1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class StrainGradientDamagePlasticityEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataStrainGradientDamagePlasticity1D();

    //! @brief copy constructor
    ConstitutiveStaticDataStrainGradientDamagePlasticity1D(ConstitutiveStaticDataStrainGradientDamagePlasticity1D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataStrainGradientDamagePlasticity1D* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataStrainGradientDamagePlasticity1D& operator= (ConstitutiveStaticDataStrainGradientDamagePlasticity1D const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataStrainGradientDamagePlasticity1D* AsStrainGradientDamagePlasticity1D();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataStrainGradientDamagePlasticity1D* AsStrainGradientDamagePlasticity1D()const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief local accumulated plastic strain 0 DP   1 Rankine
    FullMatrix<double,2,1> mKappa;

    //! @brief plastic strain
    EngineeringStrain1D mPlasticStrain;

    //! @brief previous nonlocal total strain
    EngineeringStrain1D mPrevNonlocalTotalStrain;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATASTRAINGRADIENTDAMAGEPLASTICITY1D_H
