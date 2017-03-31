#include "mechanics/loads/LoadSurfaceConstDirection3D.h"

using namespace NuTo;

LoadSurfaceConstDirection3D::LoadSurfaceConstDirection3D(StructureBase* rStructure, int rElementGroupId,
                                                         int rNodeGroupId, const Eigen::VectorXd& rLoadVector)
    : LoadSurfaceBase3D(rStructure, rElementGroupId, rNodeGroupId)
{
    mLoadVector = rLoadVector;
}


void LoadSurfaceConstDirection3D::CalculateSurfaceLoad(Eigen::Vector3d& rCoordinates, Eigen::Vector3d& rNormal,
                                                       Eigen::Vector3d& rLoadVector) const
{
    rLoadVector = mLoadVector;
}
