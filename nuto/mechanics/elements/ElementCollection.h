#pragma once
#include "nuto/mechanics/elements/ElementInterface.h"
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/elements/ElementFem.h"
#include "nuto/mechanics/elements/ElementIga.h"

namespace NuTo
{
//! @brief interface for all the cell operations, simply forwarding the corresponding element interfaces
//! @remark two benefits: a) avoids forwarding all the methods of ElementInterface for _both_ the
//                           coordinate element and the dof elements.
//                        b) allows storing elements (implementations of ElementInterface) by value
//                           in the class ElementCollectionImpl.
class ElementCollection
{
public:
    virtual ~ElementCollection() = default;
    virtual const ElementInterface& CoordinateElement() const = 0;
    virtual const ElementInterface& DofElement(DofType) const = 0;
    virtual const Shape& GetShape() const = 0;
};

//! @brief implementation of the interface ElementCollection for arbitrary element types that are derived from
//! ElementInterface
//! @tparam TElement element type stored by value
//! @remark This class stores elements by value. This is nice since it allows the copy/move operations, value
//! semantics. Additionally, the compiler is free to eliminate all the copies (whenever that is the right thing to do).
//! The access to the underlying elements is provided via const-reference. This reference is implicitly casted to the
//! base class ElementInterface.
template <typename TElement>
class ElementCollectionImpl : public ElementCollection
{
public:
    static_assert(std::is_base_of<ElementInterface, TElement>::value,
                  "TElement must be a descendant of ElementInterface");

    ElementCollectionImpl(TElement coordinateElement)
        : mCoordinateElement(coordinateElement)
        , mShape(coordinateElement.GetShape())
    {
    }

    //! @brief adds a dof element to the collection
    //! @param dofType dof type
    //! @param dofElement element to add
    void AddDofElement(DofType dofType, TElement dofElement)
    {
        mDofElements.Insert(dofType, dofElement);
        // The alternative implementation with
        //        mDofElements[dofType] = dofElement
        // will fail. The operator[] must be able default construct a new TElement, if it does not exist. It would then
        // return a reference to it and dofElement can be copied/moved into it. Our elements are not default
        // constructable (may require nodes, interpolations, ...). Thus, use Insert here, that copies/moves the entity
        // into the DofContainer without a temporary, default constructed TElement.
    }

    //! @brief Getter for CoordinateElement
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    const TElement& CoordinateElement() const override
    {
        return mCoordinateElement;
    }

    //! @brief nonconst Getter for CoordinateElement
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    TElement& CoordinateElement()
    {
        return mCoordinateElement;
    }

    //! @brief Getter for DofElements
    //! @param dofType dof type
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    const TElement& DofElement(DofType dofType) const override
    {
        return mDofElements[dofType];
    }

    //! @brief nonconst Getter for DofElements
    //! @param dofType dof type
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    TElement& DofElement(DofType dofType)
    {
        return mDofElements.At(dofType);
    }

    bool Has(DofType dof) const
    {
        return mDofElements.Has(dof);
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    TElement mCoordinateElement;
    DofContainer<TElement> mDofElements;
    const Shape& mShape;
};

// using ElementCollectionFem = ElementCollectionImpl<NuTo::DofElementFem>;

template <int TDimParameter>
using ElementCollectionIga = ElementCollectionImpl<NuTo::ElementIga<TDimParameter>>;


//! @brief implementation of the interface ElementCollection for ElementFem (uses different element types for dofs and
//! coordinates)
//! @remark This class stores elements by value. This is nice since it allows the copy/move operations, value
//! semantics. Additionally, the compiler is free to eliminate all the copies (whenever that is the right thing to do).
//! The access to the underlying elements is provided via const-reference. This reference is implicitly casted to the
//! base class ElementInterface.
class ElementCollectionFem : public ElementCollection
{
public:
    static_assert(std::is_base_of<ElementInterface, CoordinateElementFem>::value,
                  "TElement must be a descendant of ElementInterface");
    static_assert(std::is_base_of<ElementInterface, DofElementFem>::value,
                  "TElement must be a descendant of ElementInterface");

    ElementCollectionFem(CoordinateElementFem& coordinateElement)
        : mCoordinateElement(coordinateElement)
        , mShape(coordinateElement.GetShape())
    {
    }

    //! @brief adds a dof element to the collection
    //! @param dofType dof type
    //! @param dofElement element to add
    void AddDofElement(DofType dofType, DofElementFem dofElement)
    {
        mDofElements.Insert(dofType, dofElement);
        // The alternative implementation with
        //        mDofElements[dofType] = dofElement
        // will fail. The operator[] must be able default construct a new TElement, if it does not exist. It would then
        // return a reference to it and dofElement can be copied/moved into it. Our elements are not default
        // constructable (may require nodes, interpolations, ...). Thus, use Insert here, that copies/moves the entity
        // into the DofContainer without a temporary, default constructed TElement.
    }

    //! @brief Getter for CoordinateElement
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    const CoordinateElementFem& CoordinateElement() const override
    {
        return mCoordinateElement;
    }

    //! @brief nonconst Getter for CoordinateElement
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    CoordinateElementFem& CoordinateElement()
    {
        return mCoordinateElement;
    }

    //! @brief Getter for DofElements
    //! @param dofType dof type
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    const DofElementFem& DofElement(DofType dofType) const override
    {
        return mDofElements[dofType];
    }

    //! @brief nonconst Getter for DofElements
    //! @param dofType dof type
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    DofElementFem& DofElement(DofType dofType)
    {
        return mDofElements.At(dofType);
    }

    bool Has(DofType dof) const
    {
        return mDofElements.Has(dof);
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    CoordinateElementFem& mCoordinateElement;
    DofContainer<DofElementFem> mDofElements;
    const Shape& mShape;
};

} /* NuTo */
