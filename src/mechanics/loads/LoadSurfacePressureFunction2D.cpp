
#include "mechanics/loads/LoadSurfacePressureFunction2D.h"


//! @brief constructor
NuTo::LoadSurfacePressureFunction2D::LoadSurfacePressureFunction2D(int rLoadCase,
                                                                   StructureBase* rStructure,
                                                                   int rElementGroupId,
                                                                   int rNodeGroupId,
                                                                   const std::function<Eigen::Vector2d(Eigen::Vector2d)> &rLoadFunction)
: LoadSurfaceBase2D(rLoadCase,rStructure,rElementGroupId,rNodeGroupId), mLoadFunction(rLoadFunction)
{}

void NuTo::LoadSurfacePressureFunction2D::CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates,
                                                               Eigen::Vector2d& rNormal,
                                                               Eigen::Vector2d& rLoadVector)const
{
    (void)rNormal;
    (void)rLoadVector;
    rLoadVector = mLoadFunction(rCoordinates);
}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadSurfacePressureFunction2D)
#endif
