#pragma once

namespace NuTo
{
namespace Element
{

enum class eUpdateType
{
    STATICDATA = 0, //!< @ToDo[eUpdateType]: Description
    TMPSTATICDATA, //!< @ToDo[eUpdateType]: Description
    SWITCHMULTISCALE2NONLINEAR //!< move the fine scale model in a multiscale approach to the nonlinear part
};

enum class eOutput
{
    GAP_MATRIX_MORTAR, //!< gap matrix according to mortar discretization
    INTERNAL_GRADIENT, //!<
    INTERNAL_GRADIENT_ELASTIC, //!< calculates internal gradient for the case that the state variables remain unchanged
    EXTERNAL_GRADIENT, //!< TODO: calculate external forces in element
    HESSIAN_0_TIME_DERIVATIVE, //!<
    HESSIAN_0_TIME_DERIVATIVE_ELASTIC, //!<
    HESSIAN_1_TIME_DERIVATIVE, //!<
    HESSIAN_2_TIME_DERIVATIVE, //!<
    LUMPED_HESSIAN_2_TIME_DERIVATIVE, //!<
    UPDATE_STATIC_DATA,
    UPDATE_TMP_STATIC_DATA,
    IP_DATA, //!< this is primarily for plotting, give the 3D state  so for plane stress there is a z-component in the
             //!strain
    GLOBAL_ROW_DOF, //!< calculates the row dofs of the local element matrices
    GLOBAL_COLUMN_DOF //!< calculates the column dofs of the local element matrices
};

} // namespace Element
} // namespace NuTo
