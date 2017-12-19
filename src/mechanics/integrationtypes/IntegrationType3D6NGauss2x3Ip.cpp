#include <cassert>
#include "mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.h"

//! @brief constructor
NuTo::IntegrationType3D6NGauss2x3Ip::IntegrationType3D6NGauss2x3Ip()
{
}

Eigen::VectorXd NuTo::IntegrationType3D6NGauss2x3Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    /**
     *  x and y [index 0 and 1] taken from 2D3N3IP triangle
     *  z taken from 1D2N2IP truss
     */

    assert(rIpNum >= 0 && rIpNum < 6);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector3d({1. / 6., 1. / 6., -0.577350269189626});
    case 1:
        return Eigen::Vector3d({4. / 6., 1. / 6., -0.577350269189626});
    case 2:
        return Eigen::Vector3d({1. / 6., 4. / 6., -0.577350269189626});
    case 3:
        return Eigen::Vector3d({1. / 6., 1. / 6., 0.577350269189626});
    case 4:
        return Eigen::Vector3d({4. / 6., 1. / 6., 0.577350269189626});
    case 5:
        return Eigen::Vector3d({1. / 6., 4. / 6., 0.577350269189626});
    default:
        throw Exception(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D6NGauss2x3Ip::GetNumIntegrationPoints() const
{
    return 6;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D6NGauss2x3Ip::GetIntegrationPointWeight(int) const
{
    return 1 / 6.;
}
