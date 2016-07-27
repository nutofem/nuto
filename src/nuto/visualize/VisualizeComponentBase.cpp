// $Id$ 
// VisualizeComponentBase.cpp
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

#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeException.h"

int NuTo::VisualizeComponentBase::GetElementId()const
{
    throw VisualizeException(__PRETTY_FUNCTION__, "Visualization component has no ElementId.");
}

int NuTo::VisualizeComponentBase::GetIp()const
{
    throw VisualizeException(__PRETTY_FUNCTION__, "Visualization component has no Integration point.");
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponentBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponentBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponentBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponentBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponentBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::VisualizeComponentBase)
#endif // ENABLE_SERIALIZATION
