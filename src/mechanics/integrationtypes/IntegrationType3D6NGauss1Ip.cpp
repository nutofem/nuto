#include <cassert>
#include "mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.h"


Eigen::VectorXd NuTo::IntegrationType3D6NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum == 0);
    return Eigen::Vector3d({1. / 3., 1. / 3., 0});
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D6NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D6NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 1;
}
