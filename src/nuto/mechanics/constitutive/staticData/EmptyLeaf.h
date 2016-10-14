#pragma once

#include "nuto/mechanics/constitutive/staticData/Component.h"

namespace NuTo
{
namespace Constitutive
{
namespace StaticData
{
class EmptyLeaf : public Component
{
public:
    static EmptyLeaf* Create()
    {
        return new EmptyLeaf();
    }

    EmptyLeaf* Clone() const override
    {
        return new EmptyLeaf();
    }

    void ShiftToPast() override {};
    void ShiftToFuture() override {};
    void AllocateAdditionalData(int) override {};

    void WriteComponent(SerializeStreamOut& rStream) override {}
    void ReadComponent(SerializeStreamIn& rStream) override {}
};

} // namespace StaticData
} // namespace Constitutive
} // namespace NuTo
