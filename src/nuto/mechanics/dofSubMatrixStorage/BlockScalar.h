#pragma once

#include "nuto/mechanics/dofSubMatrixStorage/BlockStorageBase.h"

namespace NuTo
{


//! @author Volker Hirthammer
//! @date February 16, 2015
//! @brief ...
class BlockScalar : public BlockStorageBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    BlockScalar() {}
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    //! @param rDofStatus ... reference to DofStatus
    BlockScalar(const DofStatus& rDofStatus);

    //! @brief copy constructor
    BlockScalar(const BlockScalar&  rOther) = default;

#ifndef SWIG

    //! @brief move constructor
    //! @param rOther ... other BlockScalar
    BlockScalar(      BlockScalar&& rOther) = default;

    //! @brief destructor
    ~BlockScalar() = default;



    // Operator overloads
    // ------------------


    //! @brief copy assignment operator (only active DOFs)
    //! @param rOther ... other BlockScalar
    BlockScalar& operator=(const BlockScalar&  rOther);


    //! @brief move assignment operator
    //! @param rOther ... other BlockScalar
    BlockScalar& operator=(      BlockScalar&& rOther) = default;


    //! @brief == operator (only active DOFs)
    //! @param rOther ... other BlockScalar
    bool operator==(const NuTo::BlockScalar& rOther) const;

    //! @brief != operator (only active DOFs)
    //! @param rOther ... other BlockScalar
    bool operator!=(const NuTo::BlockScalar& rOther) const;

    //! @brief [] operator
    //! @param rDof ... degree of freedom
    //! @return reference to map entry, adds the entry to map, if needed
    double& operator[](Node::eDof rDof);

    //! @brief [] operator
    //! @param rDof ... degree of freedom
    //! @return reference to map entry
    const double operator[](Node::eDof rDof) const;


    //! @brief *= operator (only active DOFs)
    //! @param rRhs ... right operand
    BlockScalar& operator*=(double rRhs);

    //! @brief /= operator (only active DOFs)
    //! @param rRhs ... right operand
    BlockScalar& operator/=(double rRhs);

#ifndef SWIG
    //! @brief << operator, prints the BlockScalar
    //! @param rOstream ... output stream
    //! @param rRhs
    //! @return output steam containing the BlockScalar info
    friend std::ostream& operator<<(std::ostream& rOstream, const BlockScalar& rRhs);
#endif

    // Member functions
    // ----------------

    //! @brief defines a default value to all uninitialized, active dof types
    //! @param rDefaultValue ... default value
    void DefineDefaultValueToIninitializedDofTypes(double rDefaultValue);

    //! @brief Check if each active DOF value of this BlockScalar is smaller than the DOF value from the other BlockScalar
    //! @param rOther ... other BlockScalar
    //! @return true or false
    bool CheckDofWiseLessActivDofs(const BlockScalar& rOther) const;

    //! @brief Check if each active DOF value of this BlockScalar is greater than the DOF value from the other BlockScalar
    //! @param rOther ... other BlockScalar
    //! @return true or false
    bool CheckDofWiseGreaterActivDofs(const BlockScalar& rOther) const;

    //! @brief gets the number of columns of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of columns
    virtual int GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const override;

    //! @brief gets the number of rows of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of rows
    virtual int GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const override;

    //! @brief * operator (only active DOFs)
    //! @param rLhs ... left operand
    //! @param rRhs ... right operand
    //! @return multiplication result
    friend BlockScalar operator*(BlockScalar rLhs, double rRhs) {return rLhs *= rRhs;}
    friend BlockScalar operator*(double rLhs, BlockScalar rRhs) {return rRhs *= rLhs;}

    //! @brief / operator (only active DOFs)
    //! @param rLhs ... left operand
    //! @param rRhs ... right operand
    //! @return division result
    friend BlockScalar operator/(BlockScalar rLhs, double rRhs) {return rLhs /= rRhs;}

    friend bool operator <(const BlockScalar& rLhs, const BlockScalar& rRhs)
    {
        return rLhs.CheckDofWiseLessActivDofs(rRhs);
    }

    friend bool operator >(const BlockScalar& rLhs, const BlockScalar& rRhs)
    {
        return rLhs.CheckDofWiseGreaterActivDofs(rRhs);
    }

#endif

    //! @brief info method
    virtual void Info() const override;

private:


    // Class Member
    // ------------

    std::unordered_map<Node::eDof,double, Node::eDofHash> mData;
};

} //namespace NuTo

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::BlockScalar)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
