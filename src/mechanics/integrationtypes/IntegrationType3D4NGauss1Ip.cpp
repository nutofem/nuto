#include <cassert>
#include "mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"

//! @brief constructor
NuTo::IntegrationType3D4NGauss1Ip::IntegrationType3D4NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType3D4NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum == 0);
    return Eigen::Vector3d({0.25, 0.25, 0.25});
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D4NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D4NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 1 / 6.;
}
