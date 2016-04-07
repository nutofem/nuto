// $Id: ConstraintEnum.h 342 2010-10-18 12:39:08Z arnold2 $
#ifndef CONSTRAINTENUM_H_
#define CONSTRAINTENUM_H_

namespace NuTo
{
namespace Constraint
{
enum eDof
{
    COORDINATES=0,
    DISPLACEMENTS,
    ROTATIONS,
    TEMPERATURE
};

enum eSolutionProcedure
{
    GAUSSELIMINATION,      //!< eliminate the DOFs using Gauss elimination
    LAGRANGEMULTIPLIER     //!< use a lagrange multiplier for each constraint
};

enum eEquationSign
{
    EQUAL,       //!< constraint equation
    GREATER,     //!< x_1>RHS
    SMALLER     //!< x_1<RHS
};
}
}
#endif /* CONSTRAINTENUM_H_ */
