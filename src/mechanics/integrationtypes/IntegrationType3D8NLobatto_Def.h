#pragma once

#include "mechanics/integrationtypes/IntegrationType3D.h"

namespace NuTo
{
//! @author Jörg F. Unger, BAM
//! @date November 2013
//! @brief ... integration types in 3D with 8 nodes Lobatto integration with variable ip number
template <int T>
class IntegrationType3D8NLobatto : public IntegrationType3D
{

public:
    //! @brief constructor
    IntegrationType3D8NLobatto();

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

protected:
    //! @brief ... integration points coordinates
    std::vector<Eigen::Vector3d> mCoordinates;

    //! @brief ... weights for the integration
    std::vector<double> mWeights;


};
}

