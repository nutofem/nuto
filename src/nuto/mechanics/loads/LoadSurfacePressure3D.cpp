// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/loads/LoadSurfacePressure3D.h"


//! @brief constructor
NuTo::LoadSurfacePressure3D::LoadSurfacePressure3D(int rLoadCase, StructureBase* rStructure,int rElementGroupId,int rNodeGroupId,
		double rPressure) :
		LoadSurfaceBase3D(rLoadCase,rStructure,rElementGroupId,rNodeGroupId)
{
	mPressure = rPressure;
}

//! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
//! @param rCoordinates ... global coordinates
//! @param rNormal ... normal to the surface (pointing outwards)
//! @param rLoadVector ... load vector
void NuTo::LoadSurfacePressure3D::CalculateSurfaceLoad(NuTo::FullVector<double,3>& rCoordinates,NuTo::FullVector<double,3>& rNormal,
		NuTo::FullVector<double,3>& rLoadVector)const
{
	assert(std::abs(rNormal.Norm()-1.)<1e-5);
	rLoadVector = rNormal*(-mPressure);
}
