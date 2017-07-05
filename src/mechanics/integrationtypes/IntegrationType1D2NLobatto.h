#pragma once

#include <vector>
#include "mechanics/integrationtypes/IntegrationType1D.h"

namespace NuTo
{
//! @author Philipp Mueller, BAM
//! @date Jun 2017
//! @brief ... integration types in 1D with Lobatto integration
class   IntegrationType1D2NLobatto : public IntegrationType1D
{

public:

    //! @brief constructor
    IntegrationType1D2NLobatto(int numIps);

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

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const override;
#endif // ENABLE_VISUALIZE
private:
    //! @brief ... integration points coordinates
    std::vector<double> mIPts;
    //! @brief ... weights for the integration
    std::vector<double> mWeights;
};
} // namespace

