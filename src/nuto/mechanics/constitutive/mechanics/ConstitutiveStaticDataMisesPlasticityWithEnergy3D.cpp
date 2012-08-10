// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticityWithEnergy3D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::ConstitutiveStaticDataMisesPlasticityWithEnergy3D() :
ConstitutiveStaticDataMisesPlasticity3D() ,ConstitutiveStaticDataPrevEngineeringStressStrain3D()
{
}

//! @brief assignment operator
NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D& NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::operator= (ConstitutiveStaticDataMisesPlasticityWithEnergy3D const& rOther)
{
    NuTo::ConstitutiveStaticDataMisesPlasticity3D::operator= (rOther);
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::operator= (rOther);
    return (*this);
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS)
	{
		if (rElementType==NuTo::Element::BRICK8N || rElementType==NuTo::Element::TETRAHEDRON4N || rElementType==NuTo::Element::TETRAHEDRON10N)
			return true;
		else
			return false;
	}
	else
		return false;
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataMisesPlasticityWithEnergy3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataMisesPlasticity3D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain3D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataMisesPlasticityWithEnergy3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D)
#endif // ENABLE_SERIALIZATION
