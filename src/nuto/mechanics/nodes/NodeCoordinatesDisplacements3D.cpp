// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements3D.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacements3D::NodeCoordinatesDisplacements3D():
    NuTo::NodeCoordinates3D::NodeCoordinates3D (), NuTo::NodeDisplacements3D::NodeDisplacements3D()
{
#ifdef ENABLE_SERIALIZATION
//    boost::serialization::void_cast_register<NodeCoordinatesDisplacements3D, NodeBase>(0,0);
#endif
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacements3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacements3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacements3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacements3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates3D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements3D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacements3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacements3D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacements3D, track_always)
#endif // ENABLE_SERIALIZATION
