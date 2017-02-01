// $Id: ConstraintLinear.h 328 2010-10-01 14:39:32Z unger3 $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include "mechanics/constraints/ConstraintBase.h"

namespace NuTo
{
class ConstraintLagrange;
template<class T> class SparseMatrix;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all constraint equations
class ConstraintLinear : public ConstraintBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintLinear();

    //! @brief destructor
    ~ConstraintLinear();

    //! @brief true, if the constraint is linear
    bool IsLinear()const override
    {
        return true;
    }

    //! @brief cast to linear constraint
    NuTo::ConstraintLinear* AsConstraintLinear() override
    {
        return this;
    }

    //! @brief cast to linear constraint
    const NuTo::ConstraintLinear* AsConstraintLinear()const override
    {
        return this;
    }

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int GetNumLinearConstraints()const override = 0;

    using ConstraintBase::GetRHS;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    virtual void GetRHS(int& curConstraintEquation, Eigen::VectorXd& rRHS)const=0;

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    virtual void AddToConstraintMatrix(int& curConstraintEquation,
                                       NuTo::SparseMatrix<double>& rConstraintMatrix)const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinear)
#endif // ENABLE_SERIALIZATION
