#include "mechanics/integrationtypes/IntegrationTypeTet4.h"

const std::vector<Eigen::Vector4d> NuTo::IntegrationTypeTet4::quadratureData = {
        Eigen::Vector4d(-0.50000000000000, -0.500000000000000, -0.500000000000000, -0.105244444444444),
        Eigen::Vector4d(-0.85714285714285, -0.857142857142857, -0.857142857142857, 0.060977777777778),
        Eigen::Vector4d(-0.85714285714285, -0.857142857142857, 0.571428571428571, 0.060977777777778),
        Eigen::Vector4d(-0.85714285714285, 0.571428571428571, -0.857142857142857, 0.060977777777778),
        Eigen::Vector4d(0.571428571428571, -0.857142857142857, -0.857142857142857, 0.060977777777778),
        Eigen::Vector4d(-0.20119284766640, -0.201192847666402, -0.798807152333598, 0.199111111111111),
        Eigen::Vector4d(-0.20119284766640, -0.798807152333598, -0.201192847666402, 0.199111111111111),
        Eigen::Vector4d(-0.79880715233359, -0.201192847666402, -0.201192847666402, 0.199111111111111),
        Eigen::Vector4d(-0.20119284766640, -0.798807152333598, -0.798807152333598, 0.199111111111111),
        Eigen::Vector4d(-0.79880715233359, -0.201192847666402, -0.798807152333598, 0.199111111111111),
        Eigen::Vector4d(-0.79880715233359, -0.798807152333598, -0.201192847666402, 0.199111111111111)};


Eigen::VectorXd NuTo::IntegrationTypeTet4::GetLocalIntegrationPointCoordinates(int i) const
{
    // transform to standard triangle (lower left is origin)
    return (quadratureData[i].head(3) + Eigen::Vector3d(1., 1., 1.)) / 2.;
}


int NuTo::IntegrationTypeTet4::GetNumIntegrationPoints() const
{
    return quadratureData.size();
}


double NuTo::IntegrationTypeTet4::GetIntegrationPointWeight(int i) const
{
    // transform to standard triangle (lower left is origin)
    return 1. / 8 * quadratureData[i][3];
}
