#pragma once

#include "mechanics/integrationtypes/IntegrationType2D.h"

namespace NuTo
{

//! @author Thomas Titscher
//! @date August 2015
//! @brief ... 6th order integration in 2D with 12 points and visualization for every integration point via voronoi cells
class IntegrationType2D3NGauss12IpDetail : public IntegrationType2D
{

public:
    //! @brief constructor
    IntegrationType2D3NGauss12IpDetail();

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
    IpCellInfo GetVisualizationCells() const override
    {
        return mIpCellInfo;
    }
#endif // ENABLE_VISUALIZE

protected:

    std::vector<Eigen::Vector2d> mIntegrationPointCoordinates;
    std::vector<double> mIntegrationPointWeights;
#ifdef ENABLE_VISUALIZE
    IntegrationTypeBase::IpCellInfo mIpCellInfo;
#endif // ENABLE_VISUALIZE
};
}

