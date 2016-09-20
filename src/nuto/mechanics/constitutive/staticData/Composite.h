#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/staticData/Component.h"

#include <vector>

namespace NuTo
{
namespace Constitutive
{
namespace StaticData
{

//! @brief Static data for combined constitutive laws.
class Composite : public Component
{
public:
    static Composite* Create();
    Composite* Clone() const override;

    void AddComponent(Component* newComponent);

    Component& GetComponent(unsigned int index);

    unsigned int GetNumComponents() const;

    void ShiftToPast() override;

    void ShiftToFuture() override;

private:
    Composite() = default;
    boost::ptr_vector<Component> mComponents;
};

} // namspace StaticData
} // namespace Constitutive
} // namespace NuTo
