// $Id$

#ifndef CONSTRAINTNODEGROUPDISPLACEMENTS2D_H
#define CONSTRAINTNODEGROUPDISPLACEMENTS2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroup.h"

namespace NuTo
{
class NodeBase;
template <class T>
class Group;
template <class T>
class FullMatrix;
//! @author Daniel Arnold, ISM
//! @date June 2010
//! @brief ... class for all displacement constraints applied to a group of nodes in 2D
class ConstraintNodeGroupDisplacements2D : public ConstraintNodeGroup
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintNodeGroupDisplacements2D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);

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
    //! @brief prescribed displacement of the node
    double mValue;
    //! @brief direction of the applied constraint (normalized)
    double mDirection[2];
};
}//namespace NuTo
#endif //CONSTRAINTNODEGROUPDISPLACEMENTS2D_H

