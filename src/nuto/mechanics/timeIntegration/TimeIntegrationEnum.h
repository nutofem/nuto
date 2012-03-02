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

}
}
#endif /* TIMEINTEGRATIONENUM_H_ */
