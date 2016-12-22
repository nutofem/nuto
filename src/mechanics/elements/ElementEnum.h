// $Id$
#pragma once

#include <map>

namespace NuTo
{
namespace Element
{
enum class eElementType
{

    BOUNDARYGRADIENTDAMAGE1D=0,                     //!< boundary element for gradient models
    BOUNDARYMOISTURETRANSPORT1D,                    //!< boundary element for moisture transport
    CONTINUUMBOUNDARYELEMENT,                       //!< boundary of continuum element
    CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE, //!< boundary of 2D element --> 1D surface with additional node as dependency
    CONTINUUMELEMENT,                               //!< continuum element, dimension should not matter
    CONTINUUMELEMENTIGA,                            //!< continuum element iga, dimension should not matter
    ELEMENT1DINXD,                                  //!< one dimensional element in 2D or 3D
    ELEMENT1DSPRING,                                //!< one dimensional spring element
    ELEMENT2D,                                      //!< two dimensional element
    ELEMENT2DINTERFACE,                             //!< two dimensional element
    CONTACTINTERFACE,                               //!< conatact interface element
    ELEMENT3D,                                      //!< three dimensional element
    PLANE2D4N,
    PLANE2D3N,


};

const std::map<eElementType, std::string> GetElementTypeMap();
std::string ElementTypeToString(eElementType rType);
eElementType ElementTypeToEnum(std::string rType);


enum class eUpdateType
{
    STATICDATA=0,		          //!< @ToDo[eUpdateType]: Description
    TMPSTATICDATA,		          //!< @ToDo[eUpdateType]: Description
    SWITCHMULTISCALE2NONLINEAR    //!< move the fine scale model in a multiscale approach to the nonlinear part
};

enum class eOutput
{
    GAP_MATRIX_MORTAR,               //!< gap matrix according to mortar discretization
	INTERNAL_GRADIENT,               //!<
	INTERNAL_GRADIENT_ELASTIC,       //!< calculates internal gradient for the case that the state variables remain unchanged
	EXTERNAL_GRADIENT,               //!< TODO: calculate external forces in element
	HESSIAN_0_TIME_DERIVATIVE,       //!<
	HESSIAN_0_TIME_DERIVATIVE_ELASTIC,	//!<
	HESSIAN_1_TIME_DERIVATIVE,       //!<
	HESSIAN_2_TIME_DERIVATIVE,       //!<
	LUMPED_HESSIAN_2_TIME_DERIVATIVE,//!<
	UPDATE_STATIC_DATA,
	UPDATE_TMP_STATIC_DATA,
	IP_DATA,                         //!< this is primarily for plotting, give the 3D state  so for plane stress there is a z-component in the strain
	GLOBAL_ROW_DOF,                  //!< calculates the row dofs of the local element matrices
	GLOBAL_COLUMN_DOF                //!< calculates the column dofs of the local element matrices
};

}// namespace Element
}// namespace NuTo
