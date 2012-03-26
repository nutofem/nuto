// $Id$
// VisualizeComponentDisplacement.cpp
// created Apr 27, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentRotation.h"
#include "nuto/visualize/VisualizeException.h"

NuTo::VisualizeComponentRotation::VisualizeComponentRotation() : VisualizeComponentBase::VisualizeComponentBase()
{}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponentRotation::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentRotation::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentRotation::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentRotation::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentRotation::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentRotation::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponentRotation::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponentRotation" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(VisualizeComponentBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponentRotation" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponentRotation)
#endif // ENABLE_SERIALIZATION
