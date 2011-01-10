// $Id: ConstraintNodeFineScaleDisplacements2D.h 314 2010-09-27 16:31:43Z unger3 $

#ifndef ConstraintNodeFineScaleDisplacements2D_H
#define ConstraintNodeFineScaleDisplacements2D_H

#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/constraints/ConstraintNode.h"

namespace NuTo
{
template <class T>
class FullMatrix;
//! @author Jörg F. Unger, NU
//! @date October 2010
//! @brief ... abstract class for all constraints applied to a single node
class ConstraintLinearNodeFineScaleDisplacements2D : public ConstraintNode, public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLinearNodeFineScaleDisplacements2D(const NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const;

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    NuTo::ConstraintLinear* AsConstraintLinear();

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    const NuTo::ConstraintLinear* AsConstraintLinear()const;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    void SetRHS(double rRHS);

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    //! @param rRHS right hand side of the constraint equation
    void AddToConstraintMatrix(int& curConstraintEquation,
                               NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
                               NuTo::FullMatrix<double>& rRHS)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintLinearNodeFineScaleDisplacements2D(){};
    //! @brief prescribed displacement of the node
    double mValue;
    //! @brief direction of the applied constraint (normalized)
    double mDirection[2];
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearNodeFineScaleDisplacements2D)
#endif // ENABLE_SERIALIZATION

#endif //ConstraintNodeFineScaleDisplacements2D_H

