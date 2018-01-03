#include "mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"

Eigen::VectorXd NuTo::IntegrationType3D4NGauss1Ip::GetLocalIntegrationPointCoordinates(int) const
{
    return Eigen::Vector3d({0.25, 0.25, 0.25});
}


int NuTo::IntegrationType3D4NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}


double NuTo::IntegrationType3D4NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 1 / 6.;
}
