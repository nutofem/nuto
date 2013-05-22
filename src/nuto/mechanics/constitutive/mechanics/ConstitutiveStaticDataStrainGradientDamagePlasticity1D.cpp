// $Id: ConstitutiveStaticDataStrainGradientDamagePlasticity1D.cpp 612 2012-08-13 07:31:23Z unger3 $
// ConstitutiveStaticDataStrainGradientDamagePlasticity1D.cpp
// created May 6, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataStrainGradientDamagePlasticity1D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::ConstitutiveStaticDataStrainGradientDamagePlasticity1D()
   : NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::ConstitutiveStaticDataPrevEngineeringStressStrain1D()
{
    mKappa.setZero();
	mPlasticStrain[0] = 0.;
	mPrevNonlocalTotalStrain(0) = 0;
}

NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D& NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::operator= (NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D const& rOther)
{
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain1D::operator= (rOther);
    mKappa = rOther.mKappa;
    mPlasticStrain = rOther.mPlasticStrain;
    mPrevNonlocalTotalStrain = rOther.mPrevNonlocalTotalStrain;
    return *this;
}

NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::Clone()const
{
	return new ConstitutiveStaticDataStrainGradientDamagePlasticity1D(*this);
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataStrainGradientDamagePlasticity1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain1D)
       & BOOST_SERIALIZATION_NVP(mKappa)
       & BOOST_SERIALIZATION_NVP(mPlasticStrain)
       & BOOST_SERIALIZATION_NVP(mPrevNonlocalTotalStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataStrainGradientDamagePlasticity1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D)
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage1d static data
NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::AsStrainGradientDamagePlasticity1D()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D* NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::AsStrainGradientDamagePlasticity1D()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataStrainGradientDamagePlasticity1D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS)
	{
		if (rElementType==NuTo::Element::TRUSS1D2N || rElementType==NuTo::Element::TRUSS1D3N)
			return true;
		else
			return false;
	}
	else
		return false;
}

