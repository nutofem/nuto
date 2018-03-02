#pragma once

#include "mechanics/integrationtypes/IntegrationType2D.h"

namespace NuTo
{
//! @author Philip Huschke
//! @date July 2015
//! @brief fourth order integration in 2D with six points
class IntegrationType2D3NGauss6Ip : public IntegrationType2D
{

public:
    //! @brief constructor
    IntegrationType2D3NGauss6Ip();

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
};
}
