
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBondStress.h"
#include "nuto/visualize/VisualizeException.h"

NuTo::VisualizeComponentBondStress::VisualizeComponentBondStress() : VisualizeComponentBase::VisualizeComponentBase()
{}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponentBondStress::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBondStress::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBondStress::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBondStress::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBondStress::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBondStress::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponentBondStress::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponentBondStress" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(VisualizeComponentBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponentBondStress" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponentBondStress)
#endif // ENABLE_SERIALIZATION
