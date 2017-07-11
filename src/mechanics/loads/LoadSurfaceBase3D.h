#pragma once

#include "mechanics/MechanicsException.h"
#include "mechanics/loads/LoadBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include <vector>

namespace NuTo
{
class IntegrationTypeBase;
class NodeBase;
template <int TDim>
class ContinuumElement;
class StructureBase;

//! @brief Abstract class for all surface loads in 3D
class LoadSurfaceBase3D : public LoadBase
{

public:
    //! @brief constructor
    LoadSurfaceBase3D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

    //! @brief Calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    virtual void CalculateSurfaceLoad(Eigen::Vector3d& rCoordinates, Eigen::Vector3d& rNormal,
                                      Eigen::Vector3d& rLoadVector) const = 0;

protected:
    std::vector<std::pair<const ContinuumElement<3>*, int>> mVolumeElements;
    IntegrationTypeBase* mIntegrationTypeTriangleGauss1;
    IntegrationTypeBase* mIntegrationTypeTriangleGauss2;

    IntegrationTypeBase* mIntegrationTypeQuadGauss1;
    IntegrationTypeBase* mIntegrationTypeQuadGauss2;
    IntegrationTypeBase* mIntegrationTypeQuadLobatto2;
    IntegrationTypeBase* mIntegrationTypeQuadLobatto3;
    IntegrationTypeBase* mIntegrationTypeQuadLobatto4;
};
} // namespace NuTo
