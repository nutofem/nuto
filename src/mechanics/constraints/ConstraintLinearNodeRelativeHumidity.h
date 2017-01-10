#pragma once

#include "mechanics/constraints/ConstraintLinear.h"
#include "mechanics/constraints/ConstraintNode.h"

namespace NuTo
{

class ConstraintLinearNodeRelativeHumidity : public ConstraintLinear, public ConstraintNode
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
                                ConstraintLinearNodeRelativeHumidity            (const NodeBase* rNode, double rValue);


    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    virtual void                AddToConstraintMatrix                           (int& curConstraintEquation, NuTo::SparseMatrix<double>& rConstraintMatrix) const override;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int                 GetNumLinearConstraints                         () const override;

    //!@brief returns the right hand side of the constraint equations
    //!@return rRHS
    virtual double              GetRHS                                          () const override;

    //!@brief sets/modifies the right hand side of the constraint equation
    //!@param rRHS new right hand side
    virtual void                SetRHS                                          (double rRHS) override;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param rCurConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    virtual void                GetRHS                                          (int& rCurConstraintEquation, Eigen::VectorXd& rRHS) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void                Info                                            (unsigned short rVerboseLevel) const override;

    //! @brief determines the dof type affected by the constraint
    //! @return dof type
    Node::eDof GetDofType() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast std::map containing the old and new Addresses
    virtual void SetNodePtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mNodeMapCast) override;
#endif // ENABLE_SERIALIZATION


protected:
    //! @brief ... just for serialize
                                ConstraintLinearNodeRelativeHumidity()                                                 {}

    double mRHS = 0.0;

};


} //namespace NuTo




