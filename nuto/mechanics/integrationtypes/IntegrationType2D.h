#pragma once

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
class IntegrationType2D : public IntegrationTypeBase
{

public:
    //! @brief constructor
    IntegrationType2D(){};

    int GetDimension() const override
    {
        return 2;
    }
};
}
