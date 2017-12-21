#include <cassert>
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"

//! @brief constructor
NuTo::IntegrationType2D3NGauss12Ip::IntegrationType2D3NGauss12Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 12);

    const double a = 0.063089104491502;
    const double b = 0.249286745170910;
    const double c = 0.310352451033785;
    const double d = 0.053145049844816;

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
    case 6:
        return Eigen::Vector2d({c, d});
    case 7:
        return Eigen::Vector2d({d, c});
    case 8:
        return Eigen::Vector2d({1 - c - d, c});
    case 9:
        return Eigen::Vector2d({1 - c - d, d});
    case 10:
        return Eigen::Vector2d({c, 1 - c - d});
    case 11:
        return Eigen::Vector2d({d, 1 - c - d});
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss12Ip::GetNumIntegrationPoints() const
{
    return 12;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss12Ip::GetIntegrationPointWeight(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 12);
    const double e = 0.025422453185103;
    const double f = 0.058393137863189;
    const double g = 0.041425537809187;
    switch (rIpNum)
    {
    case 0:
        return e;
        break;
    case 1:
        return e;
        break;
    case 2:
        return e;
        break;
    case 3:
        return f;
        break;
    case 4:
        return f;
        break;
    case 5:
        return f;
        break;
    case 6:
        return g;
        break;
    case 7:
        return g;
        break;
    case 8:
        return g;
        break;
    case 9:
        return g;
        break;
    case 10:
        return g;
        break;
    case 11:
        return g;
        break;
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}
