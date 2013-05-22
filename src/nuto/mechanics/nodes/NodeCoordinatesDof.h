// $Id$

#ifndef NodeCoordinatesDof_H
#define NodeCoordinatesDof_H
#include "nuto/mechanics/nodes/NodeCoordinatesDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::NodeCoordinatesDof()
: NuTo::NodeCoordinates<TNumCoordinates>::NodeCoordinates(),
  NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::NodeDof()
{
}

//! @brief ... destructor
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::~NodeCoordinatesDof()
{
}

//! @brief assignment operator
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
void NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::
     operator=(NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain> const& rOther)
{
	throw MechanicsException("NuTo::NodeCoordinatesDof");
}


//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
std::string NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::
GetNodeTypeStr()const
{
	throw MechanicsException("");
    return std::string("[NuTo::NodeDof::GetNodeTypeStr]to be done");
}


#ifdef ENABLE_VISUALIZE
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
void NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::
Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{

}
#endif // ENABLE_VISUALIZE

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>*
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::Clone()const
{
    return new NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>(*this);
}

#endif //NodeCoordinatesDof_H

