// $Id$

#include <boost/assign.hpp>

#include "nuto/mechanics/nodes/MapNodeKeyNodeEnum.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::MapNodeKeyNodeEnum::MapNodeKeyNodeEnum()
{

	//KEY = mNumCoord,mNumTimeDerivative,mNumDisp,mNumRotations,mNumTemperature,mNonlocalData,mRadius,mGrid;
	//mMapNodeKeyNodeEnum = boost::assign::list_of(NodeKey(0,1,1),Node::NodeCoordinatesTimeDerivative_0_Displacements1D);
/*	        (NodeKey(0,2,2),Node::NodeCoordinatesTimeDerivative_0_Displacements2D)
	        (NodeKey(0,3,3),Node::NodeCoordinatesTimeDerivative_0_Displacements1D)
	        (NodeKey(0,3,0,0,0,1),Node::NodeCoordinatesTimeDerivative_0_Temperature3D)
	;
	*/
	//mMapNodeKeyNodeEnum[NodeKey(0,1,1)] = Node::NodeCoordinatesTimeDerivative_0_Displacements1D;
	//boost::assign::insert( mMapNodeKeyNodeEnum )
	//  (NodeKey(0,1,1),Node::NodeCoordinatesDof);
}

//! @brief destructor
NuTo::MapNodeKeyNodeEnum::~MapNodeKeyNodeEnum ()
{

}

NuTo::Node::eNodeType NuTo::MapNodeKeyNodeEnum::GetEnum(const NodeKey& rKey)const
{
	auto it = mMapNodeKeyNodeEnum.find(rKey);
	if (it==mMapNodeKeyNodeEnum.end())
		throw MechanicsException("[NuTo::MapNodeKeyNodeEnum::GetEnum] this combination of attributes for nodes is not implemented.");
	else
		return it->second;
}
