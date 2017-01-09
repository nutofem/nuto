// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $

#include "mechanics/loads/LoadSurfaceConstDirection3D.h"


//! @brief constructor
NuTo::LoadSurfaceConstDirection3D::LoadSurfaceConstDirection3D(int rLoadCase, StructureBase* rStructure,int rElementGroupId,int rNodeGroupId,
		const Eigen::VectorXd& rLoadVector) :
		LoadSurfaceBase3D(rLoadCase,rStructure,rElementGroupId,rNodeGroupId)
{
	mLoadVector = rLoadVector;
}

//! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
//! @param rCoordinates ... global coordinates
//! @param rNormal ... normal to the surface (pointing outwards)
//! @param rLoadVector ... load vector
void NuTo::LoadSurfaceConstDirection3D::CalculateSurfaceLoad(Eigen::Vector3d& rCoordinates,
															 Eigen::Vector3d& rNormal,
															 Eigen::Vector3d& rLoadVector)const
{
	rLoadVector = mLoadVector;
}
