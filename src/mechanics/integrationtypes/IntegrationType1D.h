#pragma once

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
class IntegrationType1D : public IntegrationTypeBase
{

public:
    //! @brief constructor
    IntegrationType1D()
    {
    }

    int GetDimension() const override
    {
        return 1;
    }
};
} // namespace nuto
