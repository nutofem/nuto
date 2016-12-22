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

#include "optimize/CallbackHandler.h"

// serializes the class
#ifdef ENABLE_SERIALIZATION
template void NuTo::CallbackHandler::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::CallbackHandler::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::CallbackHandler::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::CallbackHandler::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::CallbackHandler::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::CallbackHandler::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::CallbackHandler::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize CallbackHandler" << std::endl;
#endif
ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize CallbackHandler" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::CallbackHandler)
#endif // ENABLE_SERIALIZATION
