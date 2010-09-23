
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements2D.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacements2D::NodeCoordinatesDisplacements2D():
    NuTo::NodeCoordinates2D::NodeCoordinates2D (), NuTo::NodeDisplacements2D::NodeDisplacements2D()
{
#ifdef ENABLE_SERIALIZATION
    boost::serialization::void_cast_register<NodeCoordinatesDisplacements2D, NodeBase>(0,0);
#endif
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacements2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacements2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacements2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacements2D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacements2D, track_always)
#endif // ENABLE_SERIALIZATION
