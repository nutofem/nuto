// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeVelocities1D.h"

// default constructor
NuTo::NodeVelocities1D::NodeVelocities1D() : NodeBase ()
{
	this->mVelocities[0] = 0.0;
}

// constructor
NuTo::NodeVelocities1D::NodeVelocities1D(const double rVelocities[1]) : NodeBase ()
{
	this->mVelocities[0] = rVelocities[0];
}

// returns the number of velocities of the node
int NuTo::NodeVelocities1D::GetNumVelocities()const
{
	return 1;
}

// set the velocities
void NuTo::NodeVelocities1D::SetVelocities1D(const double rVelocities[1])
{
	this->mVelocities[0] = rVelocities[0];
}

// writes the velocities of a node to the prescribed pointer
void NuTo::NodeVelocities1D::GetVelocities1D(double rVelocities[1])const
{
	rVelocities[0] = this->mVelocities[0];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeVelocities1D::GetNodeTypeStr()const
{
	return std::string("NodeVelocities1D");
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeVelocities1D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mVelocities);
}
#endif // ENABLE_SERIALIZATION
