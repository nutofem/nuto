
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsVelocitiesAccelerations1D.h"

// constructor
NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::NodeCoordinatesDisplacementsVelocitiesAccelerations1D()
:NodeCoordinates1D(),NodeDisplacements1D(),NodeVelocities1D(),NodeAccelerations1D()
{
}

// sets the global dofs
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::SetGlobalDofs(int& rDOF)
{
    NodeCoordinates1D::SetGlobalDofs(rDOF);
    NodeDisplacements1D::SetGlobalDofs(rDOF);
}

// write dof values to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates1D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    NodeDisplacements1D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
}

// extract dof values from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates1D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    NodeDisplacements1D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
}

// renumber the global dofs according to predefined ordering
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    NodeCoordinates1D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    NodeDisplacements1D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
}

// write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates1D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // set velocities based on the displacement dofs
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->NodeDisplacements1D::mDOF[0];
	assert(dof >= 0);
	double value;
	if (dof >= rActiveDofValues.GetNumRows())
	{
		dof -= rActiveDofValues.GetNumRows();
		assert(dof < rDependentDofValues.GetNumRows());
		value = rDependentDofValues(dof,0);
	}
	else
	{
		assert(dof < rActiveDofValues.GetNumRows());
		value = rActiveDofValues(dof,0);
	}
	this->NodeVelocities1D::mVelocities[0] = value;
}

// extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates1D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // based on displacement dof number write velocities into global vectors
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->NodeDisplacements1D::mDOF[0];
	assert(dof >= 0);
	double value = this->NodeVelocities1D::mVelocities[0];
	if (dof >= rActiveDofValues.GetNumRows())
	{
		dof -= rActiveDofValues.GetNumRows();
		assert(dof < rDependentDofValues.GetNumRows());
		rDependentDofValues(dof,0) = value;
	}
	else
	{
		assert(dof < rActiveDofValues.GetNumRows());
		rActiveDofValues(dof,0) = value;
    }
}

// write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates1D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // set accelerations based on the displacement dofs
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->NodeDisplacements1D::mDOF[0];
	assert(dof >= 0);
	double value;
	if (dof >= rActiveDofValues.GetNumRows())
	{
		dof -= rActiveDofValues.GetNumRows();
		assert(dof < rDependentDofValues.GetNumRows());
		value = rDependentDofValues(dof,0);
	}
	else
	{
		assert(dof < rActiveDofValues.GetNumRows());
		value = rActiveDofValues(dof,0);
	}
	this->NodeAccelerations1D::mAccelerations[0] = value;
}

// extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates1D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // based on displacement dof number write accelerations into global vectors
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->NodeDisplacements1D::mDOF[0];
	assert(dof >= 0);
	double value = this->NodeAccelerations1D::mAccelerations[0];
	if (dof >= rActiveDofValues.GetNumRows())
	{
		dof -= rActiveDofValues.GetNumRows();
		assert(dof < rDependentDofValues.GetNumRows());
		rDependentDofValues(dof,0) = value;
	}
	else
	{
		assert(dof < rActiveDofValues.GetNumRows());
		rActiveDofValues(dof,0) = value;
    }
}

// returns the type of the node
std::string NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations1D::GetNodeTypeStr()const
{
	return std::string("NodeCoordinatesDisplacementsVelocitiesAccelerations1D");
}
