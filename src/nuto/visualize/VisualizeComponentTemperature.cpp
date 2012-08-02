// $Id: VisualizeComponentTemperature.cpp 342 2010-10-18 12:39:08Z arnold2 $
// VisualizeComponentTemperature.cpp
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

#include "nuto/visualize/VisualizeComponentTemperature.h"
#include "nuto/visualize/VisualizeException.h"

NuTo::VisualizeComponentTemperature::VisualizeComponentTemperature() : VisualizeComponentBase::VisualizeComponentBase()
{}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponentTemperature::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentTemperature::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentTemperature::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentTemperature::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentTemperature::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentTemperature::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponentTemperature::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponentTemperature" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(VisualizeComponentBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponentTemperature" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponentTemperature)
#endif // ENABLE_SERIALIZATION
