// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticityWithEnergy3D.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::ConstitutiveStaticDataMisesPlasticityWithEnergy3D() :
ConstitutiveStaticDataMisesPlasticity3D() ,ConstitutiveStaticDataPrevEngineeringStressStrain3D()
{
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataMisesPlasticity3D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain3D);
}
#endif // ENABLE_SERIALIZATION
