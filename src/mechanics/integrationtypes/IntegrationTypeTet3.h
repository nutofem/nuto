#pragma once

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
//! @brief ... integration types for the triangle
class IntegrationTypeTet3 : public IntegrationTypeBase
{

public:
    static const std::vector<Eigen::Vector4d> quadratureData;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int i) const override;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int i) const override;

    int GetDimension() const override
    {
        return 3;
    }
};


} // namespace NuTo
