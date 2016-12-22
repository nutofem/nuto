// $Id: ConstraintEnum.h 342 2010-10-18 12:39:08Z arnold2 $
#pragma once

namespace NuTo
{
namespace Constraint
{
enum class eDof
{
    COORDINATES=0,
    DISPLACEMENTS,
    ROTATIONS,
    TEMPERATURE
};

enum class eSolutionProcedure
{
    GAUSSELIMINATION,      //!< eliminate the DOFs using Gauss elimination
    LAGRANGEMULTIPLIER     //!< use a lagrange multiplier for each constraint
};

enum class eEquationSign
{
    EQUAL,       //!< constraint equation
    GREATER,     //!< x_1>RHS
    SMALLER     //!< x_1<RHS
};
}// namespace Constraint
}// namespace NuTo
