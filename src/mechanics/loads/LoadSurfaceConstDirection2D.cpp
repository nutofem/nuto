#include "mechanics/loads/LoadSurfaceConstDirection2D.h"

using namespace NuTo;

LoadSurfaceConstDirection2D::LoadSurfaceConstDirection2D(StructureBase* rStructure, int rElementGroupId,
                                                         int rNodeGroupId, const Eigen::VectorXd& rLoadVector)
    : LoadSurfaceBase2D(rStructure, rElementGroupId, rNodeGroupId)
{
    mLoadVector = rLoadVector;
}


void LoadSurfaceConstDirection2D::CalculateSurfaceLoad(Eigen::Vector2d&, Eigen::Vector2d&,
                                                       Eigen::Vector2d& rLoadVector) const
{
    rLoadVector = mLoadVector;
}
