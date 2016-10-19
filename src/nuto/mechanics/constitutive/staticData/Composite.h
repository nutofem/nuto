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

    //! @brief Return a new copy of the object.
    static Composite* Create();

    //! @brief Return a new copy of the object.
    Composite* Clone() const override;

    // make component uncopyable to prevent slicing; use `Clone()`
    Composite(const Composite&) = delete;

    // make component uncopyable to prevent slicing; use `Clone()`
    Composite& operator=(const Composite&) = delete;


    void AddComponent(Component* newComponent);

    Component& GetComponent(unsigned int index);

    unsigned int GetNumComponents() const;

    void ShiftToPast() override;

    void ShiftToFuture() override;

    void AllocateAdditionalData(int numAdditionalData) override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) override;

private:
    Composite() = default;
    boost::ptr_vector<Component> mComponents;

    //! @brief defines the serialization of this class
    //! @param rStream serialize input/output stream
    template <typename TStream>
    void SerializeComposite(TStream &rStream);
};

} // namspace StaticData
} // namespace Constitutive
} // namespace NuTo
