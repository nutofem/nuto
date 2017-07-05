#pragma once

#include <functional>

#include "mechanics/loads/LoadSurfaceBase2D.h"

namespace NuTo
{
class StructureBase;

//! @brief Class for surface loads in 2D with a function as load
class LoadSurfacePressureFunction2D : public LoadSurfaceBase2D
{

public:
    //! @brief Constructor
    LoadSurfacePressureFunction2D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
                                  const std::function<Eigen::Vector2d(Eigen::Vector2d)>& rLoadFunction);

    //! @brief Calculates the surface load as a function of the coordinates
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    void CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates, Eigen::Vector2d& rNormal,
                              Eigen::Vector2d& rLoadVector) const override;

protected:
    std::function<Eigen::Vector2d(Eigen::Vector2d)> mLoadFunction;

private:
    //! @brief just for serialization
    LoadSurfacePressureFunction2D() = default;
};
} // namespace NuTo

