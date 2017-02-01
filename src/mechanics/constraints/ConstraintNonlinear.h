// $Id: ConstraintLagrange.h 328 2010-10-01 14:39:32Z unger3 $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION
#include <vector>
#include "mechanics/constraints/ConstraintBase.h"

namespace NuTo
{
class ConstraintLinear;
template<class T> class SparseMatrixCSRVector2Symmetric;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all constraint equations solved with Nonlinear
class ConstraintNonlinear : public ConstraintBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintNonlinear();

    //! @brief destructor
    virtual ~ConstraintNonlinear();

    //! @brief true, if the constraint is linear
    bool IsLinear()const override
    {
        return false;
    }

    //! @brief cast to nonlinear constraint
    NuTo::ConstraintNonlinear* AsConstraintNonlinear() override
    {
        return this;
    }

    //! @brief cast to nonlinear constraint
    const NuTo::ConstraintNonlinear* AsConstraintNonlinear()const override
    {
        return this;
    }

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofs ... row and column numbers in global system
    virtual void CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
            std::vector<int>& rGlobalDofs)const=0;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    virtual void CalculateGradientInternalPotential(Eigen::VectorXd& rResult,
            std::vector<int>& rGlobalDofs)const=0;

    //! @brief calculates the internal potential
    virtual double CalculateTotalPotential()const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintNonlinear)
#endif // ENABLE_SERIALIZATION
