#pragma once

#include <memory>

namespace NuTo
{
namespace Integrands
{
class Base
{
public:
    virtual std::unique_ptr<Base> Clone() const = 0;
    virtual ~Base() = default;

    template <typename TIntegrand>
    TIntegrand& As()
    {
        return dynamic_cast<TIntegrand&>(*this);
    }

    template <typename TIntegrand>
    const TIntegrand& As() const
    {
        return dynamic_cast<const TIntegrand&>(*this);
    }
};

} /* Integrand */
} /* NuTo */
