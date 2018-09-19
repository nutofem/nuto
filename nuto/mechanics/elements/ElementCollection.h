#pragma once
#include "nuto/mechanics/elements/ElementInterface.h"
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/elements/CoordinateElementFem.h"
#include "nuto/mechanics/elements/DofElementFem.h"
#include "nuto/mechanics/elements/ElementIga.h"

namespace NuTo
{
//! @brief interface for all the cell operations, simply forwarding the corresponding element interfaces
//! @remark two benefits: a) avoids forwarding all the methods of ElementInterface for _both_ the
//!                           coordinate element and the dof elements.
//!                        b) allows storing elements (implementations of ElementInterface) by value
//!                           in the class ElementCollectionImpl.
class ElementCollection
{
public:
    virtual ~ElementCollection() = default;
    virtual const CoordinateElementInterface& CoordinateElement() const = 0;
    virtual const DofElementInterface& DofElement(DofType) const = 0;
    virtual const Shape& GetShape() const = 0;
};

//! @brief implementation of the interface ElementCollection for arbitrary element types that are derived from
//! ElementInterface
//! @tparam TCoordinateElement coordinate element type stored by value
//! @tparam TDofElement dof element type stored by value
//! @remark This class stores elements by value. This is nice since it allows the copy/move operations, value
//! semantics. Additionally, the compiler is free to eliminate all the copies (whenever that is the right thing to do).
//! The access to the underlying elements is provided via const-reference. This reference is implicitly casted to the
//! base class Coordinate- or DofElementInterface.
template <typename TCoordinateElement, typename TDofElement>
class ElementCollectionImpl : public ElementCollection
{
public:
    static_assert(std::is_base_of<CoordinateElementInterface, TCoordinateElement>::value,
                  "TCoordinateElement must be a descendant of CoordinateElementInterface");

    static_assert(std::is_base_of<DofElementInterface, TDofElement>::value,
                  "TDofElement must be a descendant of DofElementInterface");

    ElementCollectionImpl(TCoordinateElement coordinateElement)
        : mCoordinateElement(coordinateElement)
        , mShape(coordinateElement.GetShape())
    {
    }

    //! @brief adds a dof element to the collection
    //! @param dofType dof type
    //! @param dofElement element to add
    void AddDofElement(DofType dofType, TDofElement dofElement)
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
    const TCoordinateElement& CoordinateElement() const override
    {
        return mCoordinateElement;
    }

    //! @brief nonconst Getter for CoordinateElement
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    TCoordinateElement& CoordinateElement()
    {
        return mCoordinateElement;
    }

    //! @brief Getter for DofElements
    //! @param dofType dof type
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    const TDofElement& DofElement(DofType dofType) const override
    {
        return mDofElements[dofType];
    }

    //! @brief nonconst Getter for DofElements
    //! @param dofType dof type
    //! @return reference to TElement. This is implicitly casted to a reference ElementInterface when accessed via
    //! ElementCollection
    TDofElement& DofElement(DofType dofType)
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
    TCoordinateElement mCoordinateElement;
    DofContainer<TDofElement> mDofElements;
    const Shape& mShape;
};

template <int TDimParameter>
using ElementCollectionIga = ElementCollectionImpl<NuTo::ElementIga<TDimParameter>, NuTo::ElementIga<TDimParameter>>;


//! @brief implementation of the interface ElementCollection for ElementFem (uses different element types for dofs and
//! coordinates)
class ElementCollectionFem : public ElementCollection
{
public:
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
    //! @return reference to a the CoordinateElementFem. This is implicitly casted to a reference
    //! CoordinateElementInterface when accessed via ElementCollection
    const CoordinateElementFem& CoordinateElement() const override
    {
        return mCoordinateElement;
    }

    //! @brief nonconst Getter for CoordinateElement
    //! @return reference to a the CoordinateElementFem. This is implicitly casted to a reference
    //! CoordinateElementInterface when accessed via ElementCollection
    CoordinateElementFem& CoordinateElement()
    {
        return mCoordinateElement;
    }

    //! @brief Getter for DofElements
    //! @param dofType dof type
    //! @return reference to a the DofElementFem. This is implicitly casted to a reference
    //! DofElementInterface when accessed via ElementCollection
    const DofElementFem& DofElement(DofType dofType) const override
    {
        return mDofElements[dofType];
    }

    //! @brief nonconst Getter for DofElements
    //! @param dofType dof type
    //! @return reference to a the DofElementFem. This is implicitly casted to a reference
    //! DofElementInterface when accessed via ElementCollection
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
