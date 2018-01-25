#include "mechanics/integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.h"

// constructor
NuTo::IntegrationType1D2NBoundaryGauss3Ip::IntegrationType1D2NBoundaryGauss3Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{

    switch (rIpNum)
    {
    case 0:
        return Eigen::Matrix<double, 1, 1>::Constant(-1);
    case 1:
        return Eigen::Matrix<double, 1, 1>::Constant(-0.774596669241483); // -sqr(3/5)
    case 2:
        return Eigen::Matrix<double, 1, 1>::Constant(0.0);
    case 3:
        return Eigen::Matrix<double, 1, 1>::Constant(0.774596669241483); // sqr(3/5)
    default:
        throw Exception("[NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetLocalIntegrationPointCoordinates] Ip "
                        "number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetNumIntegrationPoints() const
{
    return 4;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0:
        return 0.0; // 5/9
    case 1:
        return 0.555555555555556; // 5/9
    case 2:
        return 0.888888888888889; // 8/9
    case 3:
        return 0.555555555555556; // 5/9
    default:
        throw Exception(
                "[NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}
