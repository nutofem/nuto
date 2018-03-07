#include <cassert>
#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"

NuTo::IntegrationType3D4NGauss4Ip::IntegrationType3D4NGauss4Ip()
{
    mCoordinates[0] = Eigen::Vector3d({0.13819660112501048, 0.13819660112501048, 0.13819660112501048});
    mCoordinates[1] = Eigen::Vector3d({0.58541019662496852, 0.13819660112501048, 0.13819660112501048});
    mCoordinates[2] = Eigen::Vector3d({0.13819660112501048, 0.58541019662496852, 0.13819660112501048});
    mCoordinates[3] = Eigen::Vector3d({0.13819660112501048, 0.13819660112501048, 0.58541019662496852});
}

Eigen::VectorXd NuTo::IntegrationType3D4NGauss4Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 4);
    return mCoordinates[rIpNum];
}

int NuTo::IntegrationType3D4NGauss4Ip::GetNumIntegrationPoints() const
{
    return 4;
}

double NuTo::IntegrationType3D4NGauss4Ip::GetIntegrationPointWeight(int) const
{
    return 1. / 24.;
}
