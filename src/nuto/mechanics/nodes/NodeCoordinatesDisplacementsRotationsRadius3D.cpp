// $Id:$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsRotationsRadius3D.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::NodeCoordinatesDisplacementsRotationsRadius3D():
    NuTo::NodeCoordinates3D::NodeCoordinates3D (),
    NuTo::NodeDisplacements3D::NodeDisplacements3D(),
    NuTo::NodeRotations3D::NodeRotations3D(),
    NuTo::NodeRadius::NodeRadius()
{
#ifdef ENABLE_SERIALIZATION
//    boost::serialization::void_cast_register<NodeCoordinatesDisplacementsRotations3D, NodeBase>(0,0);
#endif
}

NuTo::NodeCoordinatesDisplacementsRotationsRadius3D& NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::operator= (NuTo::NodeCoordinatesDisplacementsRotationsRadius3D const& rOther)
{
    NuTo::NodeCoordinates3D::operator= (rOther);
    NuTo::NodeDisplacements3D::operator= (rOther);
    NuTo::NodeRotations3D::operator= (rOther);
    NuTo::NodeRadius::operator= (rOther);
    return *this;
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacementsRotationsRadius3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacementsRotationsRadius3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates3D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements3D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRotations3D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRadius);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacementsRotationsRadius3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacementsRotationsRadius3D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacementsRotationsRadius3D, track_always)
#endif // ENABLE_SERIALIZATION
