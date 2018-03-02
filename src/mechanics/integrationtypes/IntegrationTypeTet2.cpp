#include "mechanics/integrationtypes/IntegrationTypeTet2.h"

const std::vector<Eigen::Vector4d> NuTo::IntegrationTypeTet2::quadratureData = {
        Eigen::Vector4d(-0.72360679774997, -0.723606797749979, -0.723606797749979, 0.333333333333333),
        Eigen::Vector4d(0.170820393249937, -0.723606797749979, -0.723606797749979, 0.333333333333333),
        Eigen::Vector4d(-0.72360679774997, 0.170820393249937, -0.723606797749979, 0.333333333333333),
        Eigen::Vector4d(-0.72360679774997, -0.723606797749979, 0.170820393249937, 0.333333333333333)};


Eigen::VectorXd NuTo::IntegrationTypeTet2::GetLocalIntegrationPointCoordinates(int i) const
{
    // transform to standard triangle (lower left is origin)
    return (quadratureData[i].head(3) + Eigen::Vector3d(1., 1., 1.)) / 2.;
}


int NuTo::IntegrationTypeTet2::GetNumIntegrationPoints() const
{
    return quadratureData.size();
}


double NuTo::IntegrationTypeTet2::GetIntegrationPointWeight(int i) const
{
    // transform to standard triangle (lower left is origin)
    return 1. / 8 * quadratureData[i][3];
}
