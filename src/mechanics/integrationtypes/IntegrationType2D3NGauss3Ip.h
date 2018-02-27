#pragma once

#include "mechanics/integrationtypes/IntegrationType2D.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date June 2010
//! @brief ... integration types in 2D with three nodes Gauss integration and 1 integration point
class IntegrationType2D3NGauss3Ip : public IntegrationType2D
{

public:
    //! @brief constructor
    IntegrationType2D3NGauss3Ip();

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
