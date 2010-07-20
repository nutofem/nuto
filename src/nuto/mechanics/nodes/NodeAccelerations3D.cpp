// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeAccelerations3D.h"

// default constructor
NuTo::NodeAccelerations3D::NodeAccelerations3D() : NodeBase ()
{
	this->mAccelerations[0] = 0.0;
	this->mAccelerations[1] = 0.0;
	this->mAccelerations[2] = 0.0;
}

// constructor
NuTo::NodeAccelerations3D::NodeAccelerations3D(const double rAccelerations[3]) : NodeBase ()
{
	this->mAccelerations[0] = rAccelerations[0];
	this->mAccelerations[1] = rAccelerations[1];
	this->mAccelerations[2] = rAccelerations[2];
}

// returns the number of accelerations of the node
int NuTo::NodeAccelerations3D::GetNumAccelerations()const
{
	return 3;
}

// set the accelerations
void NuTo::NodeAccelerations3D::SetAccelerations3D(const double rAccelerations[3])
{
	this->mAccelerations[0] = rAccelerations[0];
	this->mAccelerations[1] = rAccelerations[1];
	this->mAccelerations[2] = rAccelerations[2];
}

// writes the accelerations of a node to the prescribed pointer
void NuTo::NodeAccelerations3D::GetAccelerations3D(double rAccelerations[3])const
{
	rAccelerations[0] = this->mAccelerations[0];
	rAccelerations[1] = this->mAccelerations[1];
	rAccelerations[2] = this->mAccelerations[2];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeAccelerations3D::GetNodeTypeStr()const
{
	return std::string("NodeAccelerations3D");
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeAccelerations3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mAccelerations);
}
#endif // ENABLE_SERIALIZATION
