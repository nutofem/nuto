// $Id$ 
// VisualizeComponentNonlocalWeight.cpp
// created Apr 27, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
#include "nuto/visualize/VisualizeException.h"
#include <sstream>

NuTo::VisualizeComponentNonlocalWeight::VisualizeComponentNonlocalWeight(const ElementBase* rElement, int rElementId, int rIp) : VisualizeComponent::VisualizeComponent(NuTo::VisualizeBase::NONLOCAL_WEIGHT)
{
    mElement = rElement;
    mElementId = rElementId;
    mIp = rIp;
}

int NuTo::VisualizeComponentNonlocalWeight::GetElementId()const
{
	return mElementId;
}

const NuTo::ElementBase* NuTo::VisualizeComponentNonlocalWeight::GetElement()const
{
	return mElement;
}

int NuTo::VisualizeComponentNonlocalWeight::GetIp()const
{
	return mIp;
}

std::string NuTo::VisualizeComponentNonlocalWeight::GetComponentName()const
{
	std::stringstream out;
	out << "NonlocalWeight_Element_" << mElementId << "_Ip_" << mIp;
	return out.str();
}

NuTo::VisualizeBase::eVisualizeWhat NuTo::VisualizeComponentNonlocalWeight::GetComponentEnum(void) const
{
    return NuTo::VisualizeBase::NONLOCAL_WEIGHT;
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VisualizeComponentNonlocalWeight::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentNonlocalWeight::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentNonlocalWeight::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentNonlocalWeight::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentNonlocalWeight::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VisualizeComponentNonlocalWeight::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VisualizeComponentNonlocalWeight::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize VisualizeComponentNonlocalWeight" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(VisualizeComponent)
       & BOOST_SERIALIZATION_NVP(mElementId)
       & BOOST_SERIALIZATION_NVP(mElement)
       & BOOST_SERIALIZATION_NVP(mIp);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize VisualizeComponentNonlocalWeight" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VisualizeComponentNonlocalWeight)
#endif // ENABLE_SERIALIZATION

