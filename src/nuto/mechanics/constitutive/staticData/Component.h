#pragma once
namespace NuTo
{
class SerializeStreamIn;
class SerializeStreamOut;

namespace Constitutive
{
namespace StaticData
{

//! @brief Base class for static data.
//! 
//! The static data is an example of the composite design pattern, albeit implemented poorly. That is, there is a 
//! hierarchy for multiple constitutive laws. Each node of the tree is a component (this class), and can either be a
//! @ref Composite or a @ref Leaf. A composite is a collection of other components, a leaf is an object containing the
//! static data. Think "composites contain components, each of which could be another composite".
class Component
{
public:
    //! @brief Return a new copy of the object.
    virtual Component* Clone() const = 0;

    // make component uncopyable to prevent slicing; use `Clone()`
    Component(const Component&) = delete;
    Component& operator=(const Component&) = delete;

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    //! The current data is copied to the previous data, all others are moved.
    virtual void ShiftToPast() = 0;

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    virtual void ShiftToFuture() = 0;

    //! @brief Allocate more time steps for the static data component.
    virtual void AllocateAdditionalData(int numAdditionalData) = 0;

    virtual void WriteComponent(SerializeStreamOut& rStream) {};
    virtual void ReadComponent(SerializeStreamIn& rStream) {};

    friend NuTo::SerializeStreamOut& operator<<(NuTo::SerializeStreamOut& rStream, Component& rData);
    friend NuTo::SerializeStreamIn& operator>>(NuTo::SerializeStreamIn& rStream, Component& rData);

protected:
    //! @brief Private constructor, use @ref Create() to allocate a component on the heap.
    Component() = default;
};

} // namspace StaticData
} // namespace Constitutive
} // namespace NuTo
