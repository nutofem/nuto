// $Id$

#ifndef DAMAGE_H
#define DAMAGE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;

//! @brief ... isotropic damage variable
//! @author Jörg F. Unger, BAM
//! @date July 2012
class Damage: public ConstitutiveOutputBase
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
    Damage();

    //! @brief copy constructor
    Damage(const Damage& rDamage);

    //! @brief return output object
    Damage& GetDamage()
    {
    	return *this;
    }

    //! @brief ... get Engineering Strain
    //! @return ... Engineering Strain (exx)
    //! @sa mDeformationGradient
    const double* GetData() const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    double mDamage;
};

}

#endif // DAMAGE_H
