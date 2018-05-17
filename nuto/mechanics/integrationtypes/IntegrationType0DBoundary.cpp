#include "nuto/mechanics/integrationtypes/IntegrationType0DBoundary.h"

using namespace NuTo;

Eigen::VectorXd NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates(int) const
{
    throw Exception("[NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType0DBoundary::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0:
        return 1;
    default:
        throw Exception("[NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight] Ip number out of range.");
    }
}
