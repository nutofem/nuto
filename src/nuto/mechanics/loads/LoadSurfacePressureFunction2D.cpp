#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/loads/LoadSurfacePressureFunction2D.h"


//! @brief constructor
NuTo::LoadSurfacePressureFunction2D::LoadSurfacePressureFunction2D(int rLoadCase,
                                                                   StructureBase* rStructure,
                                                                   int rElementGroupId,
                                                                   int rNodeGroupId,
                                                                   const std::function<NuTo::FullVector<double,2>(NuTo::FullVector<double,2>)> &rLoadFunction)
: LoadSurfaceBase2D(rLoadCase,rStructure,rElementGroupId,rNodeGroupId), mLoadFunction(rLoadFunction)
{}

//! @brief calculates the surface load as a function of the coordinates
//! @param rCoordinates ... global coordinates
//! @param rNormal ... normal to the surface (pointing outwards)
//! @param rLoadVector ... load vector
void NuTo::LoadSurfacePressureFunction2D::CalculateSurfaceLoad(NuTo::FullVector<double,2>& rCoordinates,
                                                               NuTo::FullVector<double,2>& rNormal,
                                                               NuTo::FullVector<double,2>& rLoadVector)const
{
    (void)rNormal;
    (void)rLoadVector;
    rLoadVector = mLoadFunction(rCoordinates);
}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadSurfacePressureFunction2D)
#endif
