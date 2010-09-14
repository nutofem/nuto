// $Id$

#ifndef CONSTRAINTBASE_H
#define CONSTRAINTBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
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
    enum eAttributes
    {
        COORDINATES=0,
        DISPLACEMENTS,
        ROTATIONS,
        TEMPERATURES
    };
    //! @brief constructor
    ConstraintBase();

    //! @brief destructor
    virtual ~ConstraintBase();

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int GetNumConstraintEquations()const=0;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    virtual void SetRHS(double rRHS);

    //!@brief set the strain of the periodic boundary conditions
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual void SetStrain(const NuTo::FullMatrix<double>& rStrain);

    //!@brief set the strain of the periodic boundary conditions
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual void SetCrackOpening(const NuTo::FullMatrix<double>& rStrain);

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    //! @param rRHS right hand side of the constraint equation
    virtual void AddToConstraintMatrix(int& curConstraintEquation,
                                       NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
                                       NuTo::FullMatrix<double>& rRHS)const=0;

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
#endif //CONSTRAINTBASE_H

