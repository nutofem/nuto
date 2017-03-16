#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/MechanicsException.h"
#include "mechanics/sections/Section.h"

using namespace NuTo;

Section::~Section()
{
}


double Section::GetArea(double) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Section type has no cross-section area.");
}


double Section::GetThickness() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Section type has no thickness.");
}


double Section::GetCircumference() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Section type has no circumference.");
}


bool Section::IsPlaneStrain() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Section is not a plane section.");
}


void Section::Info() const
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Section::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void NuTo::Section::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void NuTo::Section::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void NuTo::Section::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void NuTo::Section::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void NuTo::Section::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void NuTo::Section::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Section" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Section" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Section)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Section)
#endif // ENABLE_SERIALIZATION
