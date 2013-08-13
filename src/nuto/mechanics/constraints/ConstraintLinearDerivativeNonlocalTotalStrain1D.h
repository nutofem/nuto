// $Id: ConstraintLinearDerivativeNonlocalTotalStrain1D.h 625 2013-04-22 16:37:11Z unger3 $

#ifndef CONSTRAINTLINEARDERIVATIVENONLOCALTOTALSTRAIN_H
#define CONSTRAINTLINEARDERIVATIVENONLOCALTOTALSTRAIN_H

#include "nuto/mechanics/constraints/ConstraintLinear.h"

namespace NuTo
{
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
                               NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix)const;

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

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief ... just for serialize
    ConstraintLinearDerivativeNonlocalTotalStrain1D(){};

    //! @brief parent element
    const ElementBase* mParentElement;

    double mLocalIpCoordinate;
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D)
#endif // ENABLE_SERIALIZATION
#endif //CONSTRAINTLINEARDERIVATIVENONLOCALTOTALSTRAIN_H

