// $Id$

#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsVelocitiesAccelerations3D.h"

// constructor
NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::NodeCoordinatesDisplacementsVelocitiesAccelerations3D()
:NodeCoordinates3D(),NodeDisplacements3D(),NodeVelocities3D(),NodeAccelerations3D()
{
}

// sets the global dofs
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::SetGlobalDofs(int& rDOF)
{
    NodeCoordinates3D::SetGlobalDofs(rDOF);
    NodeDisplacements3D::SetGlobalDofs(rDOF);
}

// write dof values to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    NodeDisplacements3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
}

// extract dof values from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    NodeDisplacements3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
}

// renumber the global dofs according to predefined ordering
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    NodeCoordinates3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    NodeDisplacements3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
}

// write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // set velocities based on the displacement dofs
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 3; dim_count++)
    {
		int dof = this->NodeDisplacements3D::mDOF[dim_count];
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
		this->NodeVelocities3D::mVelocities[dim_count] = value;
    }
}

// extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // based on displacement dof number write velocities into global vectors
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 3; dim_count++)
    {
		int dof = this->NodeDisplacements3D::mDOF[dim_count];
		assert(dof >= 0);
		double value = this->NodeVelocities3D::mVelocities[dim_count];
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
}

// write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // set accelerations based on the displacement dofs
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 3; dim_count++)
    {
		int dof = this->NodeDisplacements3D::mDOF[dim_count];
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
		this->NodeAccelerations3D::mAccelerations[dim_count] = value;
    }
}

// extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // based on displacement dof number write accelerations into global vectors
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 3; dim_count++)
    {
		int dof = this->NodeDisplacements3D::mDOF[dim_count];
		assert(dof >= 0);
		double value = this->NodeAccelerations3D::mAccelerations[dim_count];
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
}

// returns the type of the node
std::string NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations3D::GetNodeTypeStr()const
{
	return std::string("NodeCoordinatesDisplacementsVelocitiesAccelerations3D");
}
