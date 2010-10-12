// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/metamodel/Transformation.h"

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Transformation::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Transformation::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Transformation::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Transformation::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Transformation::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Transformation::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Transformation::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Transformation" << std::endl;
#endif

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Transformation" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Transformation)
#endif // ENABLE_SERIALIZATION
