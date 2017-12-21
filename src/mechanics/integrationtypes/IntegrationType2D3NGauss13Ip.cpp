#include <cassert>
#include "mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"

//! @brief constructor
NuTo::IntegrationType2D3NGauss13Ip::IntegrationType2D3NGauss13Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss13Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 13);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({0.333333333333, 0.333333333333});
    case 1:
        return Eigen::Vector2d({0.0523383720927, 0.473830813954});
    case 2:
        return Eigen::Vector2d({0.473830813954, 0.473830813954});
    case 3:
        return Eigen::Vector2d({0.473830813954, 0.0523383720927});
    case 4:
        return Eigen::Vector2d({0.655764660738, 0.172117669631});
    case 5:
        return Eigen::Vector2d({0.172117669631, 0.172117669631});
    case 6:
        return Eigen::Vector2d({0.172117669631, 0.655764660738});
    case 7:
        return Eigen::Vector2d({0.0, 0.865307354083});
    case 8:
        return Eigen::Vector2d({0.0, 0.134692645917});
    case 9:
        return Eigen::Vector2d({0.865307354083, 0.0});
    case 10:
        return Eigen::Vector2d({0.865307354083, 0.134692645917});
    case 11:
        return Eigen::Vector2d({0.134692645917, 0.0});
    case 12:
        return Eigen::Vector2d({0.134692645917, 0.865307354083});
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss13Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss13Ip::GetNumIntegrationPoints() const
{
    return 13;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss13Ip::GetIntegrationPointWeight(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 13);
    switch (rIpNum)
    {
    case 0:
        return 0.0763544833942;
        break;
    case 1:
        return 0.0490679340394;
        break;
    case 2:
        return 0.0490679340394;
        break;
    case 3:
        return 0.0490679340394;
        break;
    case 4:
        return 0.0647842146403;
        break;
    case 5:
        return 0.0647842146403;
        break;
    case 6:
        return 0.0647842146403;
        break;
    case 7:
        return 0.0136815117611;
        break;
    case 8:
        return 0.0136815117611;
        break;
    case 9:
        return 0.0136815117611;
        break;
    case 10:
        return 0.0136815117611;
        break;
    case 11:
        return 0.0136815117611;
        break;
    case 12:
        return 0.0136815117611;
        break;
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss13Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}
