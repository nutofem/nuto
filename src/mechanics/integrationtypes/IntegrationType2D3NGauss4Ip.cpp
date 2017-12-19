#include "mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"
#include <cassert>

//! @brief constructor
NuTo::IntegrationType2D3NGauss4Ip::IntegrationType2D3NGauss4Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss4Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 4);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({1. / 3., 1. / 3.});
    case 1:
        return Eigen::Vector2d({1. / 5., 1. / 5.});
    case 2:
        return Eigen::Vector2d({3. / 5., 1. / 5.});
    case 3:
        return Eigen::Vector2d({1. / 5., 3. / 5.});
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss4Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss4Ip::GetNumIntegrationPoints() const
{
    return 4;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss4Ip::GetIntegrationPointWeight(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 4);
    switch (rIpNum)
    {
    case 0:
        return -27. / 96.;
    case 1:
        return 25. / 96.;
    case 2:
        return 25. / 96.;
    case 3:
        return 25. / 96.;
    default:
        throw Exception("[NuTo::IntegrationType2D3NGauss4Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}
