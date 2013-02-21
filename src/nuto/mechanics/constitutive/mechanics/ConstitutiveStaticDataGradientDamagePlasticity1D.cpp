// $Id: ConstitutiveStaticDataGradientDamagePlasticity1D.cpp 612 2012-08-13 07:31:23Z unger3 $
// ConstitutiveStaticDataGradientDamagePlasticity1D.cpp
// created May 6, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamagePlasticity1D.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain1D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::ConstitutiveStaticDataGradientDamagePlasticity1D()
   : NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::ConstitutiveStaticDataPrevEngineeringStressStrain1D()
{
    mKappa = 0.;
    mEpsilonTotRadial = 0.;
	mEpsilonP[0] = 0.;
}

NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D& NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::operator= (NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D const& rOther)
{
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::operator= (rOther);
    mKappa = rOther.mKappa;
    mEpsilonTotRadial = rOther.mEpsilonTotRadial;
    mEpsilonP[0] = rOther.mEpsilonP[0];
    mEpsilonP[1] = rOther.mEpsilonP[1];
    return *this;
}

NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::Clone()const
{
	return new ConstitutiveStaticDataGradientDamagePlasticity1D(*this);
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataGradientDamagePlasticity1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain1D)
       & BOOST_SERIALIZATION_NVP(mKappa)
       & BOOST_SERIALIZATION_NVP(mEpsilonTotRadial)
       & BOOST_SERIALIZATION_NVP(mEpsilonP);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataGradientDamagePlasticity1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage1d static data
NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::AsGradientDamagePlasticity1D()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::AsGradientDamagePlasticity1D()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataGradientDamagePlasticity1D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS)
	{
		if (rElementType==NuTo::Element::TRUSS1D2N || rElementType==NuTo::Element::TRUSS1D3N)
			return true;
		else
			return false;
	}
	else
		return false;
}

