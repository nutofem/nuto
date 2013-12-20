// $Id$
#ifndef TIMEINTEGRATIONENUM_H_
#define TIMEINTEGRATIONENUM_H_

namespace NuTo
{
namespace TimeIntegration
{
enum eConvergenceState
{
	CONTINUE_ITERATION=0,           //!< continue iteration to converge
	CONVERGED,                      //!< converged
    DECREASE_TIME_STEP              //!< no convergence, decrease load step
};

enum eResultType
{
	TIME=0,                         //!< time
	NODE_DOF,                       //!< nodal dof (diplacement, temperature etc..)
	REACTION_DOF_GROUP              //!< sum of reaction of nodal group (reaction force, heat flow, ..
};

}
}
#endif /* TIMEINTEGRATIONENUM_H_ */
