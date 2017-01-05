// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $
#include "math/FullMatrix.h"
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
void NuTo::LoadSurfaceConstDirection3D::CalculateSurfaceLoad(NuTo::FullVector<double,3>& rCoordinates,NuTo::FullVector<double,3>& rNormal,
		NuTo::FullVector<double,3>& rLoadVector)const
{
	rLoadVector = mLoadVector;
}
