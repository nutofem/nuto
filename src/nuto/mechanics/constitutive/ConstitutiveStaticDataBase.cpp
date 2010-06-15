// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataBase::serialize(Archive & ar, const unsigned int version)
{
    std::cout << "start serialize constitutive static data base" << std::endl;
    std::cout << "finish serialize constitutive static data base" << std::endl;
}
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage2d static data
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain()
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain] Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain()const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataBase::AsNonlocalDamagePlasticity2DPlaneStrain] Static data is not of type NonlocalDamagePlasticity2DPlaneStrain.");
}
