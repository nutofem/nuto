// $Id$
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION
#include <boost/noncopyable.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionBase.h"

double NuTo::SectionBase::GetArea() const
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::GetArea] section type has no cross-section area.") ;
}

void NuTo::SectionBase::SetArea(double rArea)
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::SetArea] section type has no cross-section area.") ;
}

double NuTo::SectionBase::GetThickness() const
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::GetThickness] section type has no thickness.");
}

void NuTo::SectionBase::SetThickness(double rThickness)
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::SetThickness] section type has no thickness.");
}

void NuTo::SectionBase::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    section pointer: " << this << std::endl;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::SectionBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SectionBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::SectionBase)
#endif // ENABLE_SERIALIZATION
