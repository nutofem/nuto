// $Id$

#ifndef LOCAL_EQ_TOTAL_INELASTIC_STRAIN_H
#define LOCAL_EQ_TOTAL_INELASTIC_STRAIN_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "math/FullVector.h"
#include "mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesity;

//! @brief ... local eq total inelastic strain
//! @author Vitaliy Kindrachuk, BAM
//! @date October 2015
class LocalEqTotalInelasticStrain: public ConstitutiveOutputBase, public FullVector<double,1>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    public:

    //! @brief ... constructor
    LocalEqTotalInelasticStrain();

    //! @brief copy constructor
    LocalEqTotalInelasticStrain(const LocalEqTotalInelasticStrain& rLocalEqTotalInelasticStrain);

    //! @brief return output object
    LocalEqTotalInelasticStrain& GetLocalEqTotalInelasticStrain()
    {
    	return *this;
    }

    const double* GetData() const;

    void SetData(const double rData);


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

#endif // LOCAL_EQ_TOTAL_INELASTIC_STRAIN_H
