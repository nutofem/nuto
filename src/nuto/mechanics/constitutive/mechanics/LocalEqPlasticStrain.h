// $Id$

#ifndef LOCAL_EQ_PLASTIC_STRAIN_H
#define LOCAL_EQ_PLASTIC_STRAIN_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;

//! @brief ... nonlocal eq plastic strain (separated for rankine and drucker prager)
//! @author Jörg F. Unger, BAM
//! @date July 2012
class LocalEqPlasticStrain: public ConstitutiveOutputBase, public FullVector<double,2>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    LocalEqPlasticStrain();

    //! @brief copy constructor
    LocalEqPlasticStrain(const LocalEqPlasticStrain& rLocalEqPlasticStrain);

    //! @brief return output object
    LocalEqPlasticStrain& GetLocalEqPlasticStrain()
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

#endif // LOCAL_EQ_PLASTIC_STRAIN_H
