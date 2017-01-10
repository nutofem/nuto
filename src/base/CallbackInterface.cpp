// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "base/CallbackInterface.h"
#include "base/Exception.h"

// serializes the class
#ifdef ENABLE_SERIALIZATION
template void NuTo::CallbackInterface::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::CallbackInterface::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::CallbackInterface::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::CallbackInterface::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::CallbackInterface::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::CallbackInterface::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::CallbackInterface::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize CallbackInterface" << std::endl;
#endif

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize CallbackInterface" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::CallbackInterface)
#endif // ENABLE_SERIALIZATION

bool NuTo::CallbackInterface::Exit(NuTo::StructureBase &)
{
    throw Exception(__PRETTY_FUNCTION__, "not implemented for this callback.");
}
