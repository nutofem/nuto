#pragma once

#include "mechanics/loads/LoadSurfaceBase2D.h"

namespace NuTo
{
class NodeBase;
class StructureBase;

//! @brief Class for surface loads in 3D with a const direction and amplitude of the load
class LoadSurfacePressure2D : public LoadSurfaceBase2D
{

public:
    //! @brief Constructor
    LoadSurfacePressure2D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId, double rPressure);

    //! @brief Calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    void CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates, Eigen::Vector2d& rNormal,
                              Eigen::Vector2d& rLoadVector) const override;

protected:
    double mPressure;

private:
    //! @brief just for serialization
    LoadSurfacePressure2D() = default;
};
} // namespace NuTo

