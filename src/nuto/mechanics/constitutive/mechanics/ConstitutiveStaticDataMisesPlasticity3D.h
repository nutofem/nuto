// $Id: ConstitutiveStaticDataMisesPlasticity3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAMISESPLASTICITY3D_H
#define CONSTITUTIVESTATICDATAMISESPLASTICITY3D_H

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class ConstitutiveMisesPlasticity;
class IpDataStaticDataBase;

class ConstitutiveStaticDataMisesPlasticity3D : virtual public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveMisesPlasticity;
public:
	//! @brief constructor
	ConstitutiveStaticDataMisesPlasticity3D();

    //! @brief copy constructor
	ConstitutiveStaticDataMisesPlasticity3D(ConstitutiveStaticDataMisesPlasticity3D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    virtual ConstitutiveStaticDataMisesPlasticity3D* Clone()const
    {
    	return new ConstitutiveStaticDataMisesPlasticity3D(*this);
    }

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //! @brief assignment operator
    ConstitutiveStaticDataMisesPlasticity3D& operator= (ConstitutiveStaticDataMisesPlasticity3D const& rOther);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief accumulated plastic strain (is not always equivalent to epsilon_p)
    double mEpsilonPEq;

    //! @brief plastic strain
    double mEpsilonP[6];

    //! @brief back stress
    double mSigmaB[6];
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataMisesPlasticity3D)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAMISESPLASTICITY3D_H
