#include <cassert>
#include "mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"

//! @brief constructor
NuTo::IntegrationType2D3NGauss3Ip::IntegrationType2D3NGauss3Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss3Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 3);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({1. / 6., 1. / 6.});
    case 1:
        return Eigen::Vector2d({4. / 6., 1. / 6.});
    case 2:
        return Eigen::Vector2d({1. / 6., 4. / 6.});
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss3Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss3Ip::GetNumIntegrationPoints() const
{
    return 3;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss3Ip::GetIntegrationPointWeight(int) const
{
    return 1 / 6.;
}
