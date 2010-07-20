// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeVelocities2D.h"

// default constructor
NuTo::NodeVelocities2D::NodeVelocities2D() : NodeBase ()
{
	this->mVelocities[0] = 0.0;
	this->mVelocities[1] = 0.0;
}

// constructor
NuTo::NodeVelocities2D::NodeVelocities2D(const double rVelocities[2]) : NodeBase ()
{
	this->mVelocities[0] = rVelocities[0];
	this->mVelocities[1] = rVelocities[1];
}

// returns the number of velocities of the node
int NuTo::NodeVelocities2D::GetNumVelocities()const
{
	return 2;
}

// set the velocities
void NuTo::NodeVelocities2D::SetVelocities2D(const double rVelocities[2])
{
	this->mVelocities[0] = rVelocities[0];
	this->mVelocities[1] = rVelocities[1];
}

// writes the velocities of a node to the prescribed pointer
void NuTo::NodeVelocities2D::GetVelocities2D(double rVelocities[2])const
{
	rVelocities[0] = this->mVelocities[0];
	rVelocities[1] = this->mVelocities[1];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeVelocities2D::GetNodeTypeStr()const
{
	return std::string("NodeVelocities2D");
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeVelocities2D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mVelocities);
}
#endif // ENABLE_SERIALIZATION
