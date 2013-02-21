// $Id$

#ifndef NONLOCALDAMAGE_H
#define NONLOCALDAMAGE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;

//! @brief ... isotropic damage variable
//! @author Jörg F. Unger, BAM
//! @date July 2012
class NonlocalDamage: public ConstitutiveInputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class ConstitutiveMisesPlasticity;
    friend class NonlocalDamagePlasticity;
    friend class Multiscale;
    public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    NonlocalDamage();

    //! @brief copy constructor
    NonlocalDamage(const NonlocalDamage& rNonlocalDamage);

    //! @brief return output object
    NonlocalDamage& GetNonlocalDamage()
    {
    	return *this;
    }

    //! @brief return damage value
    double GetNonlocalDamageValue()
    {
    	return mNonlocalDamage;
    }

    //! @brief ... set damage value
    void SetNonlocalDamage(double rDamage);


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    double mNonlocalDamage;
};

}

#endif // NONLOCALDAMAGE_H
