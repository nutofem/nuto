// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $
#include "math/FullMatrix.h"
#include "mechanics/loads/LoadSurfaceConstDirection2D.h"


//! @brief constructor
NuTo::LoadSurfaceConstDirection2D::LoadSurfaceConstDirection2D(int rLoadCase, StructureBase* rStructure,int rElementGroupId,int rNodeGroupId,
		const NuTo::FullVector<double,Eigen::Dynamic>& rLoadVector) :
		LoadSurfaceBase2D(rLoadCase,rStructure,rElementGroupId,rNodeGroupId)
{
	mLoadVector = rLoadVector;
}

//! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
//! @param rCoordinates ... global coordinates
//! @param rNormal ... normal to the surface (pointing outwards)
//! @param rLoadVector ... load vector
void NuTo::LoadSurfaceConstDirection2D::CalculateSurfaceLoad(NuTo::FullVector<double,2>& rCoordinates,
                                                             NuTo::FullVector<double,2>& rNormal,
                                                             NuTo::FullVector<double,2>& rLoadVector) const
{
    rLoadVector = mLoadVector;
}
