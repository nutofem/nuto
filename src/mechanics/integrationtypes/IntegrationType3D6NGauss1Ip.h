#pragma once

#include "mechanics/integrationtypes/IntegrationType3D.h"

namespace NuTo
{
//! @author Thomas Titscher, ISM
//! @date January 2017
//! @brief ... integration types in 3D with 1 Gauss point
class IntegrationType3D6NGauss1Ip : public IntegrationType3D
{

public:
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
    IntegrationTypeBase::IpCellInfo GetVisualizationCells() const override;
#endif // ENABLE_VISUALIZE

protected:
};
}
