// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeAccelerations1D.h"

// default constructor
NuTo::NodeAccelerations1D::NodeAccelerations1D() : NodeBase ()
{
	this->mAccelerations[0] = 0.0;
}

// constructor
NuTo::NodeAccelerations1D::NodeAccelerations1D(const double rAccelerations[1]) : NodeBase ()
{
	this->mAccelerations[0] = rAccelerations[0];
}

// returns the number of accelerations of the node
int NuTo::NodeAccelerations1D::GetNumAccelerations()const
{
	return 1;
}

// set the accelerations
void NuTo::NodeAccelerations1D::SetAccelerations1D(const double rAccelerations[1])
{
	this->mAccelerations[0] = rAccelerations[0];
}

// writes the accelerations of a node to the prescribed pointer
void NuTo::NodeAccelerations1D::GetAccelerations1D(double rAccelerations[1])const
{
	rAccelerations[0] = this->mAccelerations[0];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeAccelerations1D::GetNodeTypeStr()const
{
	return std::string("NodeAccelerations1D");
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeAccelerations1D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mAccelerations);
}
#endif // ENABLE_SERIALIZATION
