#include "mechanics/loads/LoadSurfacePressure3D.h"

using namespace NuTo;

LoadSurfacePressure3D::LoadSurfacePressure3D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
                                             double rPressure)
    : LoadSurfaceBase3D(rStructure, rElementGroupId, rNodeGroupId)
{
    mPressure = rPressure;
}


void LoadSurfacePressure3D::CalculateSurfaceLoad(Eigen::Vector3d& rCoordinates, Eigen::Vector3d& rNormal,
                                                 Eigen::Vector3d& rLoadVector) const
{
    assert(std::abs(rNormal.norm() - 1.) < 1e-5);
    rLoadVector = rNormal * (-mPressure);
}
