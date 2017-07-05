#pragma once

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
class IntegrationType3D : public IntegrationTypeBase
{

public:
    //! @brief constructor
    IntegrationType3D(){};

    int GetDimension() const override
    {
        return 3;
    }
};
} // namespace nuto
