// $Id$

#ifndef NONLOCAL_EQ_STRAIN_H
#define NONLOCAL_EQ_STRAIN_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesity;

//! @brief ... nonlocal eq strain
//! @author Thomas Titscher, BAM
//! @date November 2014
class NonlocalEqStrain: public ConstitutiveInputBase, public FullVector<double,1>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientEqStrainity;
    friend class LinearElastic;
    friend class NonlocalEqStrainity;
    friend class Multiscale;
    public:

    //! @brief ... constructor
    NonlocalEqStrain();

    //! @brief copy constructor
    NonlocalEqStrain(const NonlocalEqStrain& rNonlocalEqStrain);

    //! @brief return output object
    const NonlocalEqStrain& GetNonlocalEqStrain()const
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

#endif // NONLOCAL_EQ_STRAIN_H
