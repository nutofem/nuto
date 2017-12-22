#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"

//! @brief constructor
NuTo::IntegrationType2D3NGauss1Ip::IntegrationType2D3NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum != 0)
        throw Exception(
                "[NuTo::IntegrationType2D3NGauss1Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");

    return Eigen::Vector2d({1. / 3., 1. / 3.});
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 0.5;
}
