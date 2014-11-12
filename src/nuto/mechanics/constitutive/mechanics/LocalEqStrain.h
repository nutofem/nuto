// $Id$

#ifndef LOCAL_EQ_STRAIN_H
#define LOCAL_EQ_STRAIN_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesity;

//! @brief ... local eq strain
//! @author Thomas Titscher, BAM
//! @date November 2014
class LocalEqStrain: public ConstitutiveOutputBase, public FullVector<double,1>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    public:

    //! @brief ... constructor
    LocalEqStrain();

    //! @brief copy constructor
    LocalEqStrain(const LocalEqStrain& rLocalEqStrain);

    //! @brief return output object
    LocalEqStrain& GetLocalEqStrain()
    {
    	return *this;
    }


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
};

}

#endif // LOCAL_EQ_STRAIN_H
