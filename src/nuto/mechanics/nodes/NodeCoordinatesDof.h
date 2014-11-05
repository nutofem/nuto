// $Id$

#ifndef NodeCoordinatesDof_H
#define NodeCoordinatesDof_H
#include "nuto/mechanics/nodes/NodeCoordinatesDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain>
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::NodeCoordinatesDof()
: NuTo::NodeCoordinates<TNumCoordinates>::NodeCoordinates(),
  NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::NodeDof()
{
}

//! @brief ... destructor
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain>
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::~NodeCoordinatesDof()
{
}

//! @brief assignment operator
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain>
void NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::
     operator=(NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain> const& rOther)
{
	throw MechanicsException("NuTo::NodeCoordinatesDof");
}


//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain>
std::string NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::
GetNodeTypeStr()const
{
	throw MechanicsException("");
    return std::string("[NuTo::NodeDof::GetNodeTypeStr]to be done");
}


#ifdef ENABLE_VISUALIZE
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain>
void NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::
Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{
    //add nodes to visualize
	double GlobalPointCoor[3];
	switch (this->GetNumCoordinates())
	{
	case 1:
		this->GetCoordinates1D(GlobalPointCoor);
		GlobalPointCoor[1] = 0.;
		GlobalPointCoor[2] = 0.;
		break;
	case 2:
		this->GetCoordinates2D(GlobalPointCoor);
		GlobalPointCoor[2] = 0.;
		break;
	case 3:
		this->GetCoordinates3D(GlobalPointCoor);
		break;
	default:
		throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] invalid dimension of local coordinates");
	}
	unsigned int PointId = rVisualize.AddPoint(GlobalPointCoor);

	// store data
    for (auto WhatIter = rWhat.begin();WhatIter != rWhat.end(); WhatIter++)
    {
        switch (WhatIter->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DISPLACEMENTS:
			{
				double displacements3D[3];
				switch (this->GetNumDisplacements())
				{
				case 1:
					this->GetDisplacements1D(displacements3D);
					displacements3D[1] = 0.;
					displacements3D[2] = 0.;
					break;
				case 2:
					this->GetDisplacements2D(displacements3D);
					displacements3D[2] = 0.;
					break;
				case 3:
					this->GetDisplacements3D(displacements3D);
					break;
				default:
					throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] invalid dimension of local coordinates");
				}
				rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), displacements3D);
			}
        break;
        default:
        break;
        }
    }
}
#endif // ENABLE_VISUALIZE

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures , int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain>
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>*
NuTo::NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>::Clone()const
{
    return new NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain>(*this);
}

#endif //NodeCoordinatesDof_H

