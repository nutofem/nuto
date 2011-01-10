// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeVelocities3D.h"

// default constructor
NuTo::NodeVelocities3D::NodeVelocities3D() : NodeBase ()
{
	this->mVelocities[0] = 0.0;
	this->mVelocities[1] = 0.0;
	this->mVelocities[2] = 0.0;
}

// constructor
NuTo::NodeVelocities3D::NodeVelocities3D(const double rVelocities[3]) : NodeBase ()
{
	this->mVelocities[0] = rVelocities[0];
	this->mVelocities[1] = rVelocities[1];
	this->mVelocities[2] = rVelocities[2];
}

// returns the number of velocities of the node
int NuTo::NodeVelocities3D::GetNumVelocities()const
{
	return 3;
}

// set the velocities
void NuTo::NodeVelocities3D::SetVelocities3D(const double rVelocities[3])
{
	this->mVelocities[0] = rVelocities[0];
	this->mVelocities[1] = rVelocities[1];
	this->mVelocities[2] = rVelocities[2];
}

// writes the velocities of a node to the prescribed pointer
void NuTo::NodeVelocities3D::GetVelocities3D(double rVelocities[3])const
{
	rVelocities[0] = this->mVelocities[0];
	rVelocities[1] = this->mVelocities[1];
	rVelocities[2] = this->mVelocities[2];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeVelocities3D::GetNodeTypeStr()const
{
	return std::string("NodeVelocities3D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeVelocities3D::GetNodeType()const
{
    return Node::NodeVelocities3D;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeVelocities3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mVelocities);
}
#endif // ENABLE_SERIALIZATION
