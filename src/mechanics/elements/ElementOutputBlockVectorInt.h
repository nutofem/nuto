#pragma once

#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/elements/ElementOutputBase.h"

namespace NuTo
{


//! @author Volker Hirthammer
//! @date January 25, 2016
//! @brief ...
class ElementOutputBlockVectorInt : public ElementOutputBase, public BlockFullVector<int>
{


    // Constructor / Destructor
    // ------------------------

public:
    //! @brief constructor
    //! @param rDofStatus: Reference to DofStatus needed for block matrices
    ElementOutputBlockVectorInt(const DofStatus& rDofStatus)
        : BlockFullVector<int>(rDofStatus)
    {
    }

    //! @brief copy constructor
    //! @remark ElementOutputBlockVectorInt holds heavy data, no copies allowed
    ElementOutputBlockVectorInt(const ElementOutputBlockVectorInt& rOther) = delete;

    //! @brief move constructor
    //! @param rOther ... other ElementOutputBlockVectorInt
    ElementOutputBlockVectorInt(ElementOutputBlockVectorInt&& rOther) = default;

    //! @brief destructor
    ~ElementOutputBlockVectorInt() = default;


    // Operator overloads
    // ------------------


    //! @brief copy assignment operator
    //! @remark ElementOutputBlockVectorInt holds heavy data, no copies allowed
    ElementOutputBlockVectorInt& operator=(const ElementOutputBlockVectorInt& rOther) = delete;


    //! @brief move assignment operator
    //! @param rOther ... other ElementOutputBlockVectorInt
    ElementOutputBlockVectorInt& operator=(ElementOutputBlockVectorInt&& rOther) = default;


    // Member functions
    // ----------------

    //! @brief clones the output object
    //! @return cloned output object
    virtual ElementOutputBase* Clone() const override
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] not implemented");
    }

    //! @brief Gets the output as BlockFullVector<int>
    //! @return Output as BlockFullVector<int>
    virtual BlockFullVector<int>& GetBlockFullVectorInt() override
    {
        return *this;
    }


    // Member variables
    // ----------------

public:
};


} // namespace NuTo
