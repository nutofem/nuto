#pragma once

#include "mechanics/loads/LoadSurfaceBase2D.h"

namespace NuTo
{
class NodeBase;
class StructureBase;

//! @brief Class for surface loads in 2D with a const direction and amplitude of the load
class LoadSurfaceConstDirection2D : public LoadSurfaceBase2D
{

public:
    //! @brief constructor
    LoadSurfaceConstDirection2D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
                                const Eigen::VectorXd& rLoadVector);

    //! @brief Calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    void CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates, Eigen::Vector2d& rNormal,
                              Eigen::Vector2d& rLoadVector) const override;

protected:
    Eigen::Vector2d mLoadVector;
};
} // namespace NuTo
