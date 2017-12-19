#include <cassert>
#include "mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"

//! @brief constructor
NuTo::IntegrationType2D3NGauss16Ip::IntegrationType2D3NGauss16Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss16Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 16);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({0.333333333333, 0.333333333333});
    case 1:
        return Eigen::Vector2d({0.0814148234146, 0.459292588293});
    case 2:
        return Eigen::Vector2d({0.459292588293, 0.459292588293});
    case 3:
        return Eigen::Vector2d({0.459292588293, 0.0814148234146});
    case 4:
        return Eigen::Vector2d({0.898905543366, 0.050547228317});
    case 5:
        return Eigen::Vector2d({0.050547228317, 0.050547228317});
    case 6:
        return Eigen::Vector2d({0.050547228317, 0.898905543366});
    case 7:
        return Eigen::Vector2d({0.658861384496, 0.170569307752});
    case 8:
        return Eigen::Vector2d({0.170569307752, 0.170569307752});
    case 9:
        return Eigen::Vector2d({0.170569307752, 0.658861384496});
    case 10:
        return Eigen::Vector2d({0.00839477740996, 0.728492392955});
    case 11:
        return Eigen::Vector2d({0.00839477740996, 0.263112829635});
    case 12:
        return Eigen::Vector2d({0.728492392955, 0.00839477740996});
    case 13:
        return Eigen::Vector2d({0.728492392955, 0.263112829635});
    case 14:
        return Eigen::Vector2d({0.263112829635, 0.00839477740996});
    case 15:
        return Eigen::Vector2d({0.263112829635, 0.728492392955});
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss16Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss16Ip::GetNumIntegrationPoints() const
{
    return 16;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss16Ip::GetIntegrationPointWeight(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 16);
    switch (rIpNum)
    {
    case 0:
        return 0.0721578038389;
        break;
    case 1:
        return 0.0475458171336;
        break;
    case 2:
        return 0.0475458171336;
        break;
    case 3:
        return 0.0475458171336;
        break;
    case 4:
        return 0.0162292488116;
        break;
    case 5:
        return 0.0162292488116;
        break;
    case 6:
        return 0.0162292488116;
        break;
    case 7:
        return 0.0516086852674;
        break;
    case 8:
        return 0.0516086852674;
        break;
    case 9:
        return 0.0516086852674;
        break;
    case 10:
        return 0.0136151570872;
        break;
    case 11:
        return 0.0136151570872;
        break;
    case 12:
        return 0.0136151570872;
        break;
    case 13:
        return 0.0136151570872;
        break;
    case 14:
        return 0.0136151570872;
        break;
    case 15:
        return 0.0136151570872;
        break;
    default:
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss16Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}
