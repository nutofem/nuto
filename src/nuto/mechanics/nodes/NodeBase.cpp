#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeBase.h"

//! @brief constructor
NuTo::NodeBase::NodeBase()
{}
    
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NodeBase::serialize(Archive & ar, const unsigned int version)
    {}
#endif // ENABLE_SERIALIZATION

//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeBase::GetNumCoordinates()const
{
	return 0;
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
int NuTo::NodeBase::GetNumDisplacements()const
{
	return 0;
}

//! @brief returns the number of rotations of the node
//! @return number of rotations
int NuTo::NodeBase::GetNumRotations()const
{
	return 0;
}

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
int NuTo::NodeBase::GetNumTemperatures()const
{
	return 0;
}
