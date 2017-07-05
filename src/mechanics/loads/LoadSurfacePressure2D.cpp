#include "mechanics/loads/LoadSurfacePressure2D.h"

using namespace NuTo;

LoadSurfacePressure2D::LoadSurfacePressure2D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
                                             double rPressure)
    : LoadSurfaceBase2D(rStructure, rElementGroupId, rNodeGroupId)
{
    mPressure = rPressure;
}


void LoadSurfacePressure2D::CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates, Eigen::Vector2d& rNormal,
                                                 Eigen::Vector2d& rLoadVector) const
{
    assert(std::abs(rNormal.norm() - 1.) < 1e-5);
    rLoadVector = rNormal * (-mPressure);
}


