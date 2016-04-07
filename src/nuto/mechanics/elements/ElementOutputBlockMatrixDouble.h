#pragma once

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"

namespace NuTo
{

//! @author Volker Hirthammer
//! @date January 11, 2016
//! @brief ...
class ElementOutputBlockMatrixDouble: public ElementOutputBase, public BlockFullMatrix<double>
{



    // Constructor / Destructor
    // ------------------------

public:

    //! @brief constructor
    //! @param rDofStatus: Reference to DofStatus needed for block matrices
    ElementOutputBlockMatrixDouble(const DofStatus& rDofStatus) : BlockFullMatrix<double>(rDofStatus) {}

    //! @brief copy constructor
    //! @remark ElementOutputBlockMatrixDouble holds heavy data, no copies allowed
    ElementOutputBlockMatrixDouble(const ElementOutputBlockMatrixDouble&  rOther) = delete;

    //! @brief move constructor
    //! @param rOther ... other ElementOutputBlockMatrixDouble
    ElementOutputBlockMatrixDouble(      ElementOutputBlockMatrixDouble&& rOther) = default;

    //! @brief destructor
    ~ElementOutputBlockMatrixDouble() = default;



    // Operator overloads
    // ------------------


    //! @brief copy assignment operator
    //! @remark ElementOutputBlockMatrixDouble holds heavy data, no copies allowed
    ElementOutputBlockMatrixDouble& operator=(const ElementOutputBlockMatrixDouble&  rOther) = delete;


    //! @brief move assignment operator
    //! @param rOther ... other ElementOutputBlockMatrixDouble
    ElementOutputBlockMatrixDouble& operator=(      ElementOutputBlockMatrixDouble&& rOther) = default;


    // Member functions
    // ----------------

    //! @brief clones the output object
    //! @return cloned output object
    virtual ElementOutputBase* Clone() const override
    {
        throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] not implemented");
    }

    //! @brief Gets the output as BlockFullMatrix<double>
    //! @return Output as BlockFullMatrix<double>
    virtual BlockFullMatrix<double>& GetBlockFullMatrixDouble() override
    {
        return *this;
    }



    // Member variables
    // ----------------

public:

};


} // namespace NuTo
