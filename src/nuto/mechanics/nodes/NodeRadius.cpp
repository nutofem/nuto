// $Id:$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeRadius.h"

//! @brief constructor
NuTo::NodeRadius::NodeRadius() : NodeBase ()
{
	this->mRadius=0;
}

//! @brief constructor
NuTo::NodeRadius::NodeRadius (double rRadius)  : NodeBase ()
{
	this->mRadius=rRadius;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeRadius::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeRadius::serialize(Archive & ar, const unsigned int version)
{
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
	   & BOOST_SERIALIZATION_NVP(mRadius);
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeRadius)
BOOST_CLASS_TRACKING(NuTo::NodeRadius, track_always)
#endif  // ENABLE_SERIALIZATION

//! @brief returns the number of radii
//! @return number of radii
int NuTo::NodeRadius::GetNumRadius()const
{
	return 1;
}

//! @brief returns the radius of the node
//! @param rRadius ... radius
void NuTo::NodeRadius::GetRadius(double rRadius[1])const
{
	rRadius[0]=mRadius;
}

//! @brief set the radius
//! @param rRadius  given radius
void NuTo::NodeRadius::SetRadius(const double rRadius[1])
{
	mRadius=rRadius[0];
}


