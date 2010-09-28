// $Id: ConstitutiveStaticDataNonlocalDamagePlasticity3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATANONLOCALDAMAGEPLASTICITY3D_H
#define CONSTITUTIVESTATICDATANONLOCALDAMAGEPLASTICITY3D_H

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{

class ConstitutiveStaticDataNonlocalDamagePlasticity3D : public ConstitutiveStaticDataPrevEngineeringStressStrain3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveStaticDataMisesPlasticity3D;
public:
	//! @brief constructor
    ConstitutiveStaticDataNonlocalDamagePlasticity3D();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief accumulated plastic strain scaled by nonlocal radius and a factor considering the fracture energy
    double mKappa;

    //! @brief nonlocal damage variable
    double mOmega;

    //! @brief plastic strain
    double mEpsilonP[6];
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity3D)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATANONLOCALDAMAGEPLASTICITY3D_H
