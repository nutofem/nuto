// $Id$

#ifndef CONSTRAINTNODEGROUPDISPLACEMENTS2D_H
#define CONSTRAINTNODEGROUPDISPLACEMENTS2D_H

#include "nuto/mechanics/constraints/ConstraintLinear.h"
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
class ConstraintLinearNodeGroupDisplacements2D : public ConstraintNodeGroup, public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLinearNodeGroupDisplacements2D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const;

    //!@brief sets/modifies the right hand side of the constraint equation
    //!@param rRHS new right hand side
    void SetRHS(double rRHS);

    //!@brief returns the right hand side of the constraint equation
    double GetRHS()const;

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    //! @param rRHS right hand side of the constraint equation
    void AddToConstraintMatrix(int& curConstraintEquation,
                               NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
                               NuTo::FullMatrix<double>& rRHS)const;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeGroupDisplacements2D::Info] to be implemented.");
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
    ConstraintLinearNodeGroupDisplacements2D(){};

    //! @brief direction of the applied constraint (normalized)
    double mDirection[2];

    double mRHS;
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearNodeGroupDisplacements2D)
#endif // ENABLE_SERIALIZATION
#endif //CONSTRAINTNODEGROUPDISPLACEMENTS2D_H

