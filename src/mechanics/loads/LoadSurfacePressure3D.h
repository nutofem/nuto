#pragma once

#include "mechanics/loads/LoadSurfaceBase3D.h"

namespace NuTo
{
class StructureBase;

//! @brief Class for surface loads in 3D with a const direction and amplitude of the load
class LoadSurfacePressure3D : public LoadSurfaceBase3D
{

public:
    //! @brief Constructor
    LoadSurfacePressure3D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId, double rPressure);

    //! @brief Calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    void CalculateSurfaceLoad(Eigen::Vector3d& rCoordinates, Eigen::Vector3d& rNormal,
                              Eigen::Vector3d& rLoadVector) const override;

protected:
    double mPressure;
};
} // namespace NuTo
