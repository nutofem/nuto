#ifndef CONSTRAINTLINEARNODERELATIVEHUMIDITY_H
#define CONSTRAINTLINEARNODERELATIVEHUMIDITY_H

#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/constraints/ConstraintNode.h"

namespace NuTo
{

class ConstraintLinearNodeRelativeHumidity : public ConstraintLinear, public ConstraintNode
{
public:
                                ConstraintLinearNodeRelativeHumidity            (const NodeBase* rNode, double rValue);

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    virtual void                AddToConstraintMatrix                           (int& curConstraintEquation, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix) const override;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int                 GetNumLinearConstraints                         () const override;

    //!@brief returns the right hand side of the constraint equations
    //!@return rRHS
    virtual double              GetRHS                                          () const override;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param rCurConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    virtual void                GetRHS                                          (int& rCurConstraintEquation, NuTo::FullVector<double,Eigen::Dynamic>& rRHS) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void                Info                                            (unsigned short rVerboseLevel) const override;

protected:
    double mRHS = 0.0;

};


} //namespace NuTo




#endif // CONSTRAINTLINEARNODERELATIVEHUMIDITY_H
