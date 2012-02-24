// $Id: $
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataLatticeConcrete2D.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataLatticeConcrete2D::ConstitutiveStaticDataLatticeConcrete2D() : ConstitutiveStaticDataBase()
{
	mEpsilonMax = 0.;

}

//! @brief assignment operator
NuTo::ConstitutiveStaticDataLatticeConcrete2D& NuTo::ConstitutiveStaticDataLatticeConcrete2D::operator= (ConstitutiveStaticDataLatticeConcrete2D const& rOther)
{
    NuTo::ConstitutiveStaticDataBase::operator= (rOther);
    mEpsilonMax = rOther.mEpsilonMax;

    return (*this);
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataLatticeConcrete2D::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const
{
	if (rConstitutiveType==NuTo::Constitutive::LATTICE_CONCRETE)
	{
		if (rElementType==NuTo::Element::BRICK8N || rElementType==NuTo::Element::TETRAHEDRON4N || rElementType==NuTo::Element::TETRAHEDRON10N)
			return true;
		else
			return false;
	}
	else
		return false;
}

//!@ brief reinterpret as lattice concrete 2D static data
NuTo::ConstitutiveStaticDataLatticeConcrete2D* NuTo::ConstitutiveStaticDataLatticeConcrete2D::AsConstitutiveStaticDataLatticeConcrete2D()
{
	return this;
}

//!@ brief reinterpret as lattice concrete 2D static data
const NuTo::ConstitutiveStaticDataLatticeConcrete2D* NuTo::ConstitutiveStaticDataLatticeConcrete2D::AsConstitutiveStaticDataLatticeConcrete2D()const
{
	return this;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataLatticeConcrete2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataLatticeConcrete2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mEpsilonMax);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataLatticeConcrete2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataLatticeConcrete2D)
#endif // ENABLE_SERIALIZATION
