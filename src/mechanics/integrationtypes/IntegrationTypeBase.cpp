#include <iostream>

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

void NuTo::IntegrationTypeBase::Info(int rVerboseLevel) const
{
    if (rVerboseLevel > 2)
    {
        for (int count = 0; count < GetNumIntegrationPoints(); count++)
        {
            std::cout << "    IP " << count << " weight " << GetIntegrationPointWeight(count) << std::endl;
            std::cout << "        coordinates ";
            std::cout << GetLocalIntegrationPointCoordinates(count);
        }
    }
}

Eigen::MatrixXd NuTo::IntegrationTypeBase::GetNaturalIntegrationPointCoordinates() const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented in base class.");
}
