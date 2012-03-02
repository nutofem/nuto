// $Id: ConstraintNode.h 314 2010-09-27 16:31:43Z unger3 $

#ifndef CONSTRAINTGLOBALCRACKOPENING_H
#define CONSTRAINTGLOBALCRACKOPENING_H

#include <iostream>

#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class StructureMultiscale;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all constraints applied to a single node
class ConstraintLinearGlobalCrackOpening : public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintLinearGlobalCrackOpening(const StructureMultiscale* rStructure, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
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
    void GetRHS(int& curConstraintEquation,NuTo::FullMatrix<double>& rRHS)const;

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

    //! @brief direction of the applied constraint (normalized)
    double mDirection[2];

    double mRHS;

    const StructureMultiscale* mStructure;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearGlobalCrackOpening)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTGLOBALCRACKANGLE_H

