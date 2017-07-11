// $Id $
#pragma once


#include "mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputDummy : public ElementOutputBase
{
public:
    ElementOutputDummy()
    {
    }

    ElementOutputDummy* Clone() const override
    {
        return new ElementOutputDummy(*this);
    }
};
}
