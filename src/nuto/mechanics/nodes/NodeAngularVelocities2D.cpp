// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeAngularVelocities2D.h"

// default constructor
NuTo::NodeAngularVelocities2D::NodeAngularVelocities2D() : NodeBase ()
{
	this->mAngularVelocities[0] = 0.0;
}

// constructor
NuTo::NodeAngularVelocities2D::NodeAngularVelocities2D(const double rAngularVelocities[1]) : NodeBase ()
{
	this->mAngularVelocities[0] = rAngularVelocities[0];
}

// returns the number of velocities of the node
int NuTo::NodeAngularVelocities2D::GetNumAngularVelocities()const
{
	return 1;
}

// set the velocities
void NuTo::NodeAngularVelocities2D::SetAngularVelocities2D(const double rAngularVelocities[1])
{
	this->mAngularVelocities[0] = rAngularVelocities[0];
}

// writes the velocities of a node to the prescribed pointer
void NuTo::NodeAngularVelocities2D::GetAngularVelocities2D(double rAngularVelocities[1])const
{
	rAngularVelocities[0] = this->mAngularVelocities[0];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeAngularVelocities2D::GetNodeTypeStr()const
{
	return std::string("NodeAngularVelocities2D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeAngularVelocities2D::GetNodeType()const
{
    return Node::NodeAngularVelocities2D;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeAngularVelocities2D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mAngularVelocities);
}
#endif // ENABLE_SERIALIZATION
