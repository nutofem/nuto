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

#include "nuto/mechanics/constitutive/mechanics/RelativeInterfaceDisplacements3D.h"

NuTo::RelativeInterfaceDisplacements3D::RelativeInterfaceDisplacements3D()
{
    this->mRelativeInterfaceDisplacements[0] = 0.0;
    this->mRelativeInterfaceDisplacements[1] = 0.0;
    this->mRelativeInterfaceDisplacements[2] = 0.0;
}

NuTo::RelativeInterfaceDisplacements3D::RelativeInterfaceDisplacements3D(const RelativeInterfaceDisplacements3D& rOther)
{
	mRelativeInterfaceDisplacements[0] = rOther.mRelativeInterfaceDisplacements[0];
	mRelativeInterfaceDisplacements[1] = rOther.mRelativeInterfaceDisplacements[1];
	mRelativeInterfaceDisplacements[2] = rOther.mRelativeInterfaceDisplacements[2];
}

unsigned int NuTo::RelativeInterfaceDisplacements3D::GetNumberOfComponents() const
{
    return 3;
}

const double* NuTo::RelativeInterfaceDisplacements3D::GetRelativeInterfaceDisplacements3D() const
{
    return this->mRelativeInterfaceDisplacements;
}

void NuTo::RelativeInterfaceDisplacements3D::GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements3D& rRelativeInterfaceDisplacements) const
{
	rRelativeInterfaceDisplacements = NuTo::RelativeInterfaceDisplacements3D(*this);
}

void NuTo::RelativeInterfaceDisplacements3D::SetRelativeInterfaceDisplacements3D(const double* rRelativeInterfaceDisplacements)
{
    this->mRelativeInterfaceDisplacements[0] = rRelativeInterfaceDisplacements[0];
    this->mRelativeInterfaceDisplacements[1] = rRelativeInterfaceDisplacements[1];
    this->mRelativeInterfaceDisplacements[2] = rRelativeInterfaceDisplacements[2];
}

void NuTo::RelativeInterfaceDisplacements3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    relative interface displacement components: " << this->mRelativeInterfaceDisplacements[0] << ", "
              << this->mRelativeInterfaceDisplacements[1] << ", "
              << this->mRelativeInterfaceDisplacements[2] << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::RelativeInterfaceDisplacements3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RelativeInterfaceDisplacements3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize RelativeInterfaceDisplacements3D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mRelativeInterfaceDisplacements);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize RelativeInterfaceDisplacements3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
