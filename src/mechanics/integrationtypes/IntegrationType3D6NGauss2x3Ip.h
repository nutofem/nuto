#pragma once

#include "mechanics/integrationtypes/IntegrationType3D.h"

namespace NuTo
{
//! @author Thomas Titscher, ISM
//! @date January 2017
//! @brief ... integration types in 3D with 6 nodes Gauss integration and 2x3 integration points
class IntegrationType3D6NGauss2x3Ip : public IntegrationType3D
{

public:
    //! @brief constructor
    IntegrationType3D6NGauss2x3Ip();

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

    IntegrationTypeBase::IpCellInfo GetVisualizationCells() const override;
};
}
