// $Id$

#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/mechanics/RelativeInterfaceDisplacements2D.h"

NuTo::RelativeInterfaceDisplacements2D::RelativeInterfaceDisplacements2D()
{
    this->mRelativeInterfaceDisplacements[0] = 0.0;
    this->mRelativeInterfaceDisplacements[1] = 0.0;
}

NuTo::RelativeInterfaceDisplacements2D::RelativeInterfaceDisplacements2D(const RelativeInterfaceDisplacements2D& rOther)
{
	mRelativeInterfaceDisplacements[0] = rOther.mRelativeInterfaceDisplacements[0];
	mRelativeInterfaceDisplacements[1] = rOther.mRelativeInterfaceDisplacements[1];
}

unsigned int NuTo::RelativeInterfaceDisplacements2D::GetNumberOfComponents() const
{
    return 2;
}

const double* NuTo::RelativeInterfaceDisplacements2D::GetRelativeInterfaceDisplacements2D() const
{
    return this->mRelativeInterfaceDisplacements;
}

void NuTo::RelativeInterfaceDisplacements2D::GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements2D& rRelativeInterfaceDisplacements) const
{
	rRelativeInterfaceDisplacements = NuTo::RelativeInterfaceDisplacements2D(*this);
}

void NuTo::RelativeInterfaceDisplacements2D::SetRelativeInterfaceDisplacements2D(const double* rRelativeInterfaceDisplacements)
{
    this->mRelativeInterfaceDisplacements[0] = rRelativeInterfaceDisplacements[0];
    this->mRelativeInterfaceDisplacements[1] = rRelativeInterfaceDisplacements[1];
}

void NuTo::RelativeInterfaceDisplacements2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    relative interface displacement components: " << this->mRelativeInterfaceDisplacements[0] << ", "
              << this->mRelativeInterfaceDisplacements[1] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::RelativeInterfaceDisplacements2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RelativeInterfaceDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize RelativeInterfaceDisplacements2D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mRelativeInterfaceDisplacements);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize RelativeInterfaceDisplacements2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
