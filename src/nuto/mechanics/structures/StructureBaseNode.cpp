// $Id: $

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"

//! @brief sets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeSetDisplacements(int rNode, const FullMatrix<double>& rDisplacements)
{
	NodeBase* nodePtr=NodeGetNodePtr(rNode);

	if (rDisplacements.GetNumColumns()!=1)
	throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] Displacement matrix has to have a single column.");
	try
	{
		switch (rDisplacements.GetNumRows())
		{
		case 1:
			nodePtr->SetDisplacements1D(rDisplacements.mEigenMatrix.data());
		break;
		case 2:
			nodePtr->SetDisplacements2D(rDisplacements.mEigenMatrix.data());
		break;
		case 3:
			nodePtr->SetDisplacements3D(rDisplacements.mEigenMatrix.data());
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] The number of displacement components is either 1, 2 or 3.");
		}
	}
    catch(NuTo::MechanicsException & b)
	{
    	b.AddMessage("[NuTo::StructureBase::NodeSetDisplacements] Error setting displacements.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] Error setting displacements of node (unspecified exception).");
	}
}

//! @brief gets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeGetDisplacements(int rNode, FullMatrix<double>& rDisplacements)const
{
	const NodeBase* nodePtr = NodeGetNodePtr(rNode);

	try
	{
		switch (nodePtr->GetNumDisplacements())
		{
		case 1:
			rDisplacements.Resize(1,1);
			nodePtr->GetDisplacements1D(rDisplacements.mEigenMatrix.data());
		break;
		case 2:
			rDisplacements.Resize(2,1);
			nodePtr->GetDisplacements2D(rDisplacements.mEigenMatrix.data());
		break;
		case 3:
			rDisplacements.Resize(3,1);
			nodePtr->GetDisplacements3D(rDisplacements.mEigenMatrix.data());
		break;
		case 0:
			throw MechanicsException("[NuTo::StructureBase::NodeGetDisplacements] Node has no displacements.");
		break;
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGetDisplacements] Error getting displacements.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGetDisplacements] Error getting displacements of node (unspecified exception).");
	}
}

//! @brief returns the number of global dofs
//! @return number of global dofs
int NuTo::StructureBase::NodeGetNumberGlobalDofs()const
{
    if (mNodeNumberingRequired)
        throw MechanicsException("[NuTo::StructureBase::NodeGetNumberGlobalDofs] Number the DOF's first.");
    return mNumDofs;
}

//! @brief returns the number of active dofs
//! @return number of active dofs
int NuTo::StructureBase::NodeGetNumberActiveDofs()const
{
    if (mNodeNumberingRequired)
        throw MechanicsException("[NuTo::StructureBase::NodeGetNumberActiveDofs] Number the DOF's first.");
    return mNumActiveDofs;
}

