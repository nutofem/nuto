// $Id $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentHeatFlux.h"
#include "nuto/visualize/VisualizeException.h"

NuTo::VisualizeComponentHeatFlux::VisualizeComponentHeatFlux() : VisualizeComponentBase::VisualizeComponentBase()
{}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponentHeatFlux::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentHeatFlux::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentHeatFlux::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentHeatFlux::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentHeatFlux::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentHeatFlux::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponentHeatFlux::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponentHeatFlux" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(VisualizeComponentBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponentHeatFlux" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponentHeatFlux)
#endif // ENABLE_SERIALIZATION

