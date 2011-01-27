// $Id: ConstraintNode.h 314 2010-09-27 16:31:43Z unger3 $

#ifndef CONSTRAINTGLOBALCRACKOPENING_H
#define CONSTRAINTGLOBALCRACKOPENING_H

#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class StructureIp;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all constraints applied to a single node
class ConstraintLinearGlobalCrackOpening : public ConstraintBase, public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintLinearGlobalCrackOpening(const StructureIp* rStructure, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const;

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    NuTo::ConstraintLinear* AsConstraintLinear();

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    const NuTo::ConstraintLinear* AsConstraintLinear()const;

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
        std::cout << "ConstraintLinearGlobalCrackOpening" << std::endl;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintLinearGlobalCrackOpening(){};

    //! @brief prescribed crack opening
    double mValue;
    //! @brief direction of the applied constraint (normalized)
    double mDirection[2];

    const StructureIp* mStructure;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearGlobalCrackOpening)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTGLOBALCRACKANGLE_H

