#include "nuto/mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.h"


Eigen::VectorXd NuTo::IntegrationType3D6NGauss1Ip::GetLocalIntegrationPointCoordinates(int) const
{
    return Eigen::Vector3d({1. / 3., 1. / 3., 0});
}


int NuTo::IntegrationType3D6NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}


double NuTo::IntegrationType3D6NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 1;
}
