// $Id:$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsRadius2D.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacementsRadius2D::NodeCoordinatesDisplacementsRadius2D():
    NuTo::NodeCoordinates2D::NodeCoordinates2D (), NuTo::NodeDisplacements2D::NodeDisplacements2D(), NuTo::NodeRadius::NodeRadius()
{
}

NuTo::NodeCoordinatesDisplacementsRadius2D& NuTo::NodeCoordinatesDisplacementsRadius2D::operator= (NuTo::NodeCoordinatesDisplacementsRadius2D const& rOther)
{
    NuTo::NodeCoordinates2D::operator= (rOther);
    NuTo::NodeDisplacements2D::operator= (rOther);
    NuTo::NodeRadius::operator= (rOther);
    return *this;
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacementsRadius2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacementsRadius2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRadius);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacementsRadius2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacementsRadius2D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacementsRadius2D, track_always)
#endif // ENABLE_SERIALIZATION

