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

#include "nuto/mechanics/constitutive/mechanics/RelativeInterfaceDisplacements1D.h"

NuTo::RelativeInterfaceDisplacements1D::RelativeInterfaceDisplacements1D()
{
    this->mRelativeInterfaceDisplacements = 0.0;
}

NuTo::RelativeInterfaceDisplacements1D::RelativeInterfaceDisplacements1D(const RelativeInterfaceDisplacements1D& rOther)
{
	mRelativeInterfaceDisplacements = rOther.mRelativeInterfaceDisplacements;
}

unsigned int NuTo::RelativeInterfaceDisplacements1D::GetNumberOfComponents() const
{
    return 1;
}

const double* NuTo::RelativeInterfaceDisplacements1D::GetRelativeInterfaceDisplacements1D() const
{
    return &mRelativeInterfaceDisplacements;
}

void NuTo::RelativeInterfaceDisplacements1D::GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements1D& rRelativeInterfaceDisplacements) const
{
    rRelativeInterfaceDisplacements = NuTo::RelativeInterfaceDisplacements1D(*this);
}

void NuTo::RelativeInterfaceDisplacements1D::SetRelativeInterfaceDisplacements1D(const double* rRelativeInterfaceDisplacements)
{
    this->mRelativeInterfaceDisplacements = rRelativeInterfaceDisplacements[0];
}

void NuTo::RelativeInterfaceDisplacements1D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    relative interface displacement components: " << this->mRelativeInterfaceDisplacements << std::endl;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::RelativeInterfaceDisplacements1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RelativeInterfaceDisplacements1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RelativeInterfaceDisplacements1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize RelativeInterfaceDisplacements1D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mRelativeInterfaceDisplacements);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize RelativeInterfaceDisplacements1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
