#include <cassert>
#include "mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"

//! @brief constructor
NuTo::IntegrationType2D3NGauss6Ip::IntegrationType2D3NGauss6Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss6Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 6);

    const double a = 0.445948490915965;
    const double b = 0.091576213509771;

    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({a, a});
    case 1:
        return Eigen::Vector2d({1 - 2 * a, a});
    case 2:
        return Eigen::Vector2d({a, 1 - 2 * a});
    case 3:
        return Eigen::Vector2d({b, b});
    case 4:
        return Eigen::Vector2d({1 - 2 * b, b});
    case 5:
        return Eigen::Vector2d({b, 1 - 2 * b});
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss6Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss6Ip::GetNumIntegrationPoints() const
{
    return 6;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss6Ip::GetIntegrationPointWeight(int rIpNum) const
{
    const double c = 0.111690794839005;
    const double d = 0.054975871827661;

    assert(rIpNum >= 0 && rIpNum < 6);
    switch (rIpNum)
    {
    case 0:
        return c;
    case 1:
        return c;
    case 2:
        return c;
    case 3:
        return d;
    case 4:
        return d;
    case 5:
        return d;
    default:
        throw Exception("[NuTo::IntegrationType2D3NGauss6Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}
