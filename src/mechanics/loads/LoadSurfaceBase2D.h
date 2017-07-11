#pragma once

#include "base/Exception.h"
#include "mechanics/loads/LoadBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include <vector>

namespace NuTo
{
class IntegrationTypeBase;
template <int TDim>
class ContinuumElement;
class StructureBase;

//! @brief Abstract class for all surface loads in 2D
class LoadSurfaceBase2D : public LoadBase
{

public:
    //! @brief constructor
    LoadSurfaceBase2D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

    //! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    virtual void CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates, Eigen::Vector2d& rNormal,
                                      Eigen::Vector2d& rLoadVector) const = 0;

protected:
    LoadSurfaceBase2D()
    {
    }

    std::vector<std::pair<const ContinuumElement<2>*, int>> mElements2D;
    IntegrationTypeBase* mIntegrationType2NPtr;
    IntegrationTypeBase* mIntegrationType3NPtr;
    IntegrationTypeBase* mIntegrationType4NPtr;
    IntegrationTypeBase* mIntegrationType5NPtr;
    IntegrationTypeBase* mIntegrationType3NPtrLobatto;
    IntegrationTypeBase* mIntegrationType4NPtrLobatto;
    IntegrationTypeBase* mIntegrationType5NPtrLobatto;
};
} // namespace NuTo
