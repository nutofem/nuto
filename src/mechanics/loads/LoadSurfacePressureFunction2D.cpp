#include "mechanics/loads/LoadSurfacePressureFunction2D.h"

using namespace NuTo;

LoadSurfacePressureFunction2D::LoadSurfacePressureFunction2D(
        StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
        const std::function<Eigen::Vector2d(Eigen::Vector2d)>& rLoadFunction)
    : LoadSurfaceBase2D(rStructure, rElementGroupId, rNodeGroupId)
    , mLoadFunction(rLoadFunction)
{
}


void LoadSurfacePressureFunction2D::CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates, Eigen::Vector2d& rNormal,
                                                         Eigen::Vector2d& rLoadVector) const
{
    rLoadVector = mLoadFunction(rCoordinates);
}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(LoadSurfacePressureFunction2D)
#endif
