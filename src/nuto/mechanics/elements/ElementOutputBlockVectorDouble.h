#pragma once

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{


//! @author Volker Hirthammer
//! @date January 25, 2016
//! @brief ...
class ElementOutputBlockVectorDouble: public ElementOutputBase, public BlockFullVector<double>
{



    // Constructor / Destructor
    // ------------------------

public:

    //! @brief constructor
    //! @param rDofStatus: Reference to DofStatus needed for block matrices
    ElementOutputBlockVectorDouble(const DofStatus& rDofStatus) : BlockFullVector<double>(rDofStatus) {}

    //! @brief copy constructor
    //! @remark ElementOutputBlockVectorDouble holds heavy data, no copies allowed
    ElementOutputBlockVectorDouble(const ElementOutputBlockVectorDouble&  rOther) = delete;

    //! @brief move constructor
    //! @param rOther ... other ElementOutputBlockVectorDouble
    ElementOutputBlockVectorDouble(      ElementOutputBlockVectorDouble&& rOther) = default;

    //! @brief destructor
    ~ElementOutputBlockVectorDouble() = default;



    // Operator overloads
    // ------------------


    //! @brief copy assignment operator
    //! @remark ElementOutputBlockVectorDouble holds heavy data, no copies allowed
    ElementOutputBlockVectorDouble& operator=(const ElementOutputBlockVectorDouble&  rOther) = delete;


    //! @brief move assignment operator
    //! @param rOther ... other ElementOutputBlockVectorDouble
    ElementOutputBlockVectorDouble& operator=(      ElementOutputBlockVectorDouble&& rOther) = default;


    // Member functions
    // ----------------

    //! @brief clones the output object
    //! @return cloned output object
    virtual ElementOutputBase* Clone() const override
    {
        throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] not implemented");
    }


    //! @brief Gets the output as BlockFullVector<double>
    //! @return Output as BlockFullVector<double>
    virtual BlockFullVector<double>& GetBlockFullVectorDouble() override
    {
        return *this;
    }



    // Member variables
    // ----------------

public:


};


} // namespace NuTo
