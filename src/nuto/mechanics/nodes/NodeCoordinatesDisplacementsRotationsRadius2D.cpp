// $Id:$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsRotationsRadius2D.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::NodeCoordinatesDisplacementsRotationsRadius2D():
    NuTo::NodeCoordinates2D::NodeCoordinates2D (),
    NuTo::NodeDisplacements2D::NodeDisplacements2D(),
    NuTo::NodeRotations2D::NodeRotations2D(),
    NuTo::NodeRadius::NodeRadius()
{
}

NuTo::NodeCoordinatesDisplacementsRotationsRadius2D& NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::operator= (NuTo::NodeCoordinatesDisplacementsRotationsRadius2D const& rOther)
{
    NuTo::NodeCoordinates2D::operator= (rOther);
    NuTo::NodeDisplacements2D::operator= (rOther);
    NuTo::NodeRotations2D::operator= (rOther);
    NuTo::NodeRadius::operator= (rOther);
    return *this;
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacementsRotationsRadius2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacementsRotationsRadius2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRotations2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRadius);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacementsRotationsRadius2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacementsRotationsRadius2D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacementsRotationsRadius2D, track_always)
#endif // ENABLE_SERIALIZATION

