#pragma once

#include "mechanics/integrationtypes/IntegrationType3D.h"
#include <array>

namespace NuTo
{
//! @brief ... integration types in 3D with 4 nodes Gauss integration and 1 integration points
class IntegrationType3D4NGauss4Ip : public IntegrationType3D
{

public:
    //! @brief constructor
    IntegrationType3D4NGauss4Ip();

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const override;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum) const override;

    void GetVisualizationCells(unsigned int& NumVisualizationPoints,
                               std::vector<double>& VisualizationPointLocalCoordinates,
                               unsigned int& NumVisualizationCells,
                               std::vector<NuTo::eCellTypes>& VisualizationCellType,
                               std::vector<unsigned int>& VisualizationCellsIncidence,
                               std::vector<unsigned int>& VisualizationCellsIP) const override;
protected:
    std::array<Eigen::Vector3d, 4> mCoordinates;
};
} // namespace NuTo
