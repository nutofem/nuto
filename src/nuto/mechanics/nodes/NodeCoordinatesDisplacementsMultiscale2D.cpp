// $Id: NodeCoordinatesDisplacementsMultiscale2D.cpp 343 2010-10-19 07:43:10Z arnold2 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsMultiscale2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacementsMultiscale2D::NodeCoordinatesDisplacementsMultiscale2D(NuTo::StructureMultiscale* rStructureMultiscale, bool rCrackedDomain):
    NuTo::NodeCoordinates2D::NodeCoordinates2D (), NuTo::NodeDisplacementsMultiscale2D::NodeDisplacementsMultiscale2D(rStructureMultiscale, rCrackedDomain)
{
}

NuTo::NodeCoordinatesDisplacementsMultiscale2D& NuTo::NodeCoordinatesDisplacementsMultiscale2D::operator= (NuTo::NodeCoordinatesDisplacementsMultiscale2D const& rOther)
{
    NuTo::NodeCoordinates2D::operator= (rOther);
    NuTo::NodeDisplacementsMultiscale2D::operator= (rOther);
    return *this;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacementsMultiscale2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacementsMultiscale2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacementsMultiscale2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacementsMultiscale2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacementsMultiscale2D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacementsMultiscale2D, track_always)
#endif // ENABLE_SERIALIZATION
