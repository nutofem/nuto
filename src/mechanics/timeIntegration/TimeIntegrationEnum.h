#pragma once

namespace NuTo
{
// enum class eTimeIntegrationConvergenceState
//{
//	CONTINUE_ITERATION=0,           //!< continue iteration to converge
//	CONVERGED,                      //!< converged
//    DECREASE_TIME_STEP              //!< no convergence, decrease load step
//};


enum class eTimeIntegrationResultType
{
    TIME = 0, //!< time
    NODE_DISPLACEMENT, //!< nodal diplacement
    NODE_ACCELERATION, //!< nodal accelerations
    GROUP_NODE_FORCE, //!< nodal forces
    ELEMENT_IP_STRESS, //!> ip (stress)
    ELEMENT_IP_STRAIN, //!> ip (strain)
    ELEMENT_IP_DAMAGE, //!> ip (damage)
    ELEMENT_IP_BOND_STRESS, //!> ip (bond stress)
    ELEMENT_IP_SLIP, //!> ip (slip)

};

} // namespace NuTo
