// $Id: ConstraintLinearDerivativeNonlocalTotalStrain1D.h 625 2013-04-22 16:37:11Z unger3 $

#pragma once

#include "nuto/mechanics/constraints/ConstraintLinear.h"

namespace NuTo
{
class ElementBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for all constraints with displacements applied to a single node in 1D
class ConstraintLinearDerivativeNonlocalTotalStrain1D : public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rNode ... node pointer
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLinearDerivativeNonlocalTotalStrain1D(const ElementBase* rParentElement, double rLocalIpCoordinate);

    //! @brief returns the number of constraint equations
    //! @return number of conrstraints
    int GetNumLinearConstraints()const;

    //!@brief sets/modifies the right hand side of the constraint equation
    //!@param rRHS new right hand side
    void SetRHS(double rRHS);

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void AddToConstraintMatrix(int& curConstraintEquation,
                               NuTo::SparseMatrix<double>& rConstraintMatrix)const;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const
    {
        throw MechanicsException("[NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::Info] to be implemented.");
    }

    //! @brief determines the dof type affected by the constraint
    //! @return dof type
    Node::eDof GetDofType() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes (saves) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    //! @brief deserializes (loads) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    //! @brief ElementBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the ElementBase-Pointer is done by searching and casting back the address in the map
    //! @param mNodeMapCast std::map containing the old and new addresses
    virtual void SetElementPtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mElementMapCast) override;
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief ... just for serialize
    ConstraintLinearDerivativeNonlocalTotalStrain1D(){}

    //! @brief parent element
    const ElementBase* mParentElement;

    double mLocalIpCoordinate;
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D)
#endif // ENABLE_SERIALIZATION
