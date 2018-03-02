#include "mechanics/integrationtypes/IntegrationTypeTet1.h"

const std::vector<Eigen::Vector4d> NuTo::IntegrationTypeTet1::quadratureData = {
        Eigen::Vector4d(-0.500000000000000, -0.500000000000000, -0.500000000000000, 1.333333333333333)};

Eigen::VectorXd NuTo::IntegrationTypeTet1::GetLocalIntegrationPointCoordinates(int i) const
{
    // transform to standard triangle (lower left is origin)
    return (quadratureData[i].head(3) + Eigen::Vector3d(1., 1., 1.)) / 2.;
}


int NuTo::IntegrationTypeTet1::GetNumIntegrationPoints() const
{
    return quadratureData.size();
}


double NuTo::IntegrationTypeTet1::GetIntegrationPointWeight(int i) const
{
    // transform to standard triangle (lower left is origin)
    return 1. / 8 * quadratureData[i][3];
}
