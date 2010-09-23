
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsNonlocalData2D.h"

//! @brief constructor
NuTo::NodeCoordinatesDisplacementsNonlocalData2D::NodeCoordinatesDisplacementsNonlocalData2D():
       NuTo::NodeCoordinates2D::NodeCoordinates2D (), NuTo::NodeDisplacements2D::NodeDisplacements2D(), NuTo::NodeNonlocalDataBase::NodeNonlocalDataBase()
{
    std::cout<<"call NodeCoordinatesDisplacementsNonlocalData2D constructor" << std::endl;
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinatesDisplacementsNonlocalData2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinatesDisplacementsNonlocalData2D" << std::endl;
#endif
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
//       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D)
//       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeNonlocalDataBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinatesDisplacementsNonlocalData2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinatesDisplacementsNonlocalData2D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinatesDisplacementsNonlocalData2D, track_always)
#endif // ENABLE_SERIALIZATION
