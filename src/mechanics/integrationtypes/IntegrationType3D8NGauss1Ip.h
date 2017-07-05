#pragma once

#include "mechanics/integrationtypes/IntegrationType3D.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date November 2009
//! @brief ... integration types in 3D with 8 nodes Gauss integration and 2x2x2 integration points
class IntegrationType3D8NGauss1Ip : public IntegrationType3D
{

public:
    //! @brief constructor
    IntegrationType3D8NGauss1Ip();

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


};
}

