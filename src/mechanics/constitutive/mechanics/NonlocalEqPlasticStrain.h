// $Id$

#ifndef NONLOCAL_EQ_PLASTIC_STRAIN_H
#define NONLOCAL_EQ_PLASTIC_STRAIN_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION


#include "mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;

//! @brief ... nonlocal eq plastic strain (separated for rankine and drucker prager)
//! @author Jörg F. Unger, BAM
//! @date July 2012
class NonlocalEqPlasticStrain: public ConstitutiveInputBase, public FullVector<double,2>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class GradientEqPlasticStrainPlasticity;
    friend class LinearElastic;
    friend class NonlocalEqPlasticStrainPlasticity;
    friend class Multiscale;
    public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    NonlocalEqPlasticStrain();

    //! @brief copy constructor
    NonlocalEqPlasticStrain(const NonlocalEqPlasticStrain& rNonlocalEqPlasticStrain);

    //! @brief return output object
    const NonlocalEqPlasticStrain& GetNonlocalEqPlasticStrain()const
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

#endif // NONLOCAL_EQ_PLASTIC_STRAIN_H
