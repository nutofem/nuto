#pragma once

#include "nuto/mechanics/integrationtypes/IntegrationType1D.h"

namespace NuTo
{
//! @author Stefan Eckardt, IFF
//! @date July 2010
//! @brief integration types in 1D with two nodes Gauss integration and 3 integration points
class IntegrationType1D2NBoundaryGauss3Ip : public IntegrationType1D
{

public:
    //! @brief constructor
    IntegrationType1D2NBoundaryGauss3Ip();

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
