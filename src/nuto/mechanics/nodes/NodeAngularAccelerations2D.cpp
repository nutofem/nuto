// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeAngularAccelerations2D.h"

// default constructor
NuTo::NodeAngularAccelerations2D::NodeAngularAccelerations2D() : NodeBase ()
{
	this->mAngularAccelerations[0] = 0.0;
}

// constructor
NuTo::NodeAngularAccelerations2D::NodeAngularAccelerations2D(const double rAngularAccelerations[1]) : NodeBase ()
{
	this->mAngularAccelerations[0] = rAngularAccelerations[0];
}

// returns the number of accelerations of the node
int NuTo::NodeAngularAccelerations2D::GetNumAngularAccelerations()const
{
	return 1;
}

// set the accelerations
void NuTo::NodeAngularAccelerations2D::SetAngularAccelerations2D(const double rAngularAccelerations[2])
{
	this->mAngularAccelerations[0] = rAngularAccelerations[0];
}

// writes the accelerations of a node to the prescribed pointer
void NuTo::NodeAngularAccelerations2D::GetAngularAccelerations2D(double rAngularAccelerations[2])const
{
	rAngularAccelerations[0] = this->mAngularAccelerations[0];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeAngularAccelerations2D::GetNodeTypeStr()const
{
	return std::string("NodeAngularAccelerations2D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeAngularAccelerations2D::GetNodeType()const
{
    return NuTo::Node::NodeAngularAccelerations2D;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeAngularAccelerations2D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mAngularAccelerations);
}
#endif // ENABLE_SERIALIZATION
