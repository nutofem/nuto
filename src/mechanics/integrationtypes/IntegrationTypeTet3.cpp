#include "mechanics/integrationtypes/IntegrationTypeTet3.h"

const std::vector<Eigen::Vector4d> NuTo::IntegrationTypeTet3::quadratureData = {
        Eigen::Vector4d(-0.50000000000000, -0.500000000000000, -0.500000000000000, -1.066666666666667),
        Eigen::Vector4d(-0.66666666666666, -0.666666666666667, -0.666666666666667, 0.600000000000000),
        Eigen::Vector4d(-0.66666666666666, -0.666666666666667, 0.000000000000000, 0.600000000000000),
        Eigen::Vector4d(-0.66666666666666, 0.000000000000000, -0.666666666666667, 0.600000000000000),
        Eigen::Vector4d(0.000000000000000, -0.666666666666667, -0.666666666666667, 0.600000000000000)};


Eigen::VectorXd NuTo::IntegrationTypeTet3::GetLocalIntegrationPointCoordinates(int i) const
{
    // transform to standard triangle (lower left is origin)
    return (quadratureData[i].head(3) + Eigen::Vector3d(1., 1., 1.)) / 2.;
}


int NuTo::IntegrationTypeTet3::GetNumIntegrationPoints() const
{
    return quadratureData.size();
}


double NuTo::IntegrationTypeTet3::GetIntegrationPointWeight(int i) const
{
    // transform to standard triangle (lower left is origin)
    return 1. / 8 * quadratureData[i][3];
}
