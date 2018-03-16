#pragma once

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
//! @author Thomas Titscher
//! @date Jan 2015
//! @brief integration types in 0D, more like a dummy integration type
class IntegrationType0DBoundary : public IntegrationTypeBase
{

public:
    //! @brief constructor
    IntegrationType0DBoundary() = default;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const override;

    int GetDimension() const override;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum) const override;
};
}
