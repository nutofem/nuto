// $Id$

#ifndef CONSTRAINTBASE_H
#define CONSTRAINTBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constraints/ConstraintEnum.h"

namespace NuTo
{
class EngineeringStrain2D;
class NodeBase;
class ConstraintLinear;
class ConstraintNonlinear;
class ConstraintLagrange;
template<class T> class FullMatrix;
template<class T> class SparseMatrixCSRGeneral;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all constraint equations
class ConstraintBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintBase();

    //! @brief destructor
    virtual ~ConstraintBase();

    //! @brief true, if the constraint is linear
    virtual bool IsLinear()const=0;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int GetNumLinearConstraints()const;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int GetNumLagrangeMultipliers()const;

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    virtual ConstraintLinear* AsConstraintLinear();

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    virtual const ConstraintLinear* AsConstraintLinear()const;

    //! @brief cast to nonlinear constraint - the corresponding dofs are eliminated in the global system
    virtual ConstraintNonlinear* AsConstraintNonlinear();

    //! @brief cast to nonlinear constraint - the corresponding dofs are eliminated in the global system
    virtual const ConstraintNonlinear* AsConstraintNonlinear()const;

    //! @brief cast to linear constraint - Lagrange multipliers are added to the system of equations
    virtual ConstraintLagrange* AsConstraintLagrange();

    //! @brief cast to linear constraint - Lagrange multipliers are added to the system of equations
    virtual const ConstraintLagrange* AsConstraintLagrange()const;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    virtual void SetRHS(double rRHS);

    //!@brief returns the right hand side of the constraint equations
    //!@return rRHS
    virtual double GetRHS()const;

    //!@brief set the strain of the periodic boundary conditions
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual void SetStrain(const EngineeringStrain2D& rStrain);

    //!@brief get the strain of a constrain equation
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual const EngineeringStrain2D& GetStrain()const;

    //!@brief set the strain of the periodic boundary conditions
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual void SetCrackOpening(const NuTo::FullMatrix<double>& rStrain);

    //! @brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    virtual void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr){};

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const=0;

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
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintBase)
#endif // ENABLE_SERIALIZATION
#endif //CONSTRAINTBASE_H

