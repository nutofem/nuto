// $Id: ConstitutiveStaticDataMisesPlasticity3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAMISESPLASTICITY3D_H
#define CONSTITUTIVESTATICDATAMISESPLASTICITY3D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class ConstitutiveMisesPlasticity;

class ConstitutiveStaticDataMisesPlasticity3D : virtual public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveMisesPlasticity;
public:
	//! @brief constructor
	ConstitutiveStaticDataMisesPlasticity3D();

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

#endif // CONSTITUTIVESTATICDATAMISESPLASTICITY3D_H
