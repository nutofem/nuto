// $Id$
#ifndef ELEMENTENUM_H_
#define ELEMENTENUM_H_

#include <map>
#include <algorithm>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
namespace Element
{
enum eElementType
{

    BOUNDARYGRADIENTDAMAGE1D=0,                     //!< boundary element for gradient models
    BOUNDARYMOISTURETRANSPORT1D,                    //!< boundary element for moisture transport
    CONTINUUMBOUNDARYELEMENT,                       //!< boundary of continuum element
    CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE, //!< boundary of 2D element --> 1D surface with additional node as dependency
    CONTINUUMELEMENT,                               //!< continuum element, dimension should not matter
    ELEMENT1DINXD,                                  //!< one dimensional element in 2D or 3D
    ELEMENT1DSPRING,                                //!< one dimensional spring element
    ELEMENT2D,                                      //!< two dimensional element
    ELEMENT2DINTERFACE,                             //!< two dimensional element
    INTERFACEIGA,                                   //!< interface element isogeometric analysis (e.g. contact)
    ELEMENT3D,                                      //!< three dimensional element
    PLANE2D4N,
    PLANE2D3N,

};


static inline std::map<eElementType, std::string> GetElementTypeMap()
{
    std::map<eElementType, std::string> map;
    map[CONTINUUMELEMENT]                                = "CONTINUUMELEMENT";
    map[CONTINUUMBOUNDARYELEMENT]                        = "CONTINUUMBOUNDARYELEMENT";
    map[CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE]  = "CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE";
    map[ELEMENT1DINXD]                                   = "ELEMENT1DINXD";
    map[ELEMENT2DINTERFACE]                              = "ELEMENTINTERFACE";
    map[INTERFACEIGA]                                    = "INTERFACEIGA";
    return map;
}


static inline std::string ElementTypeToString(eElementType rType)
{
    try
    {
        return GetElementTypeMap().find(rType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}

static inline eElementType ElementTypeToEnum(std::string rType)
{
    std::transform(rType.begin(), rType.end(),rType.begin(), ::toupper);

    for(auto entry : GetElementTypeMap())
        if (entry.second == rType)
            return entry.first;

    throw MechanicsException(__PRETTY_FUNCTION__, "ElementType " + rType + " has no enum equivalent or is not implemented.");
}

enum eUpdateType
{
    STATICDATA=0,		          //!< @ToDo[eUpdateType]: Description
    TMPSTATICDATA,		          //!< @ToDo[eUpdateType]: Description
    CRACK,				          //!< update for crack informations
    SWITCHMULTISCALE2NONLINEAR    //!< move the fine scale model in a multiscale approach to the nonlinear part
};

enum eOutput
{
        CONTACT_FORCE,                   //!< contact force after mortar discretization
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

}
}
#endif /* ELEMENTENUM_H_ */
