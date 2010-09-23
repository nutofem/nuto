// $Id$
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsVelocitiesAccelerations2D.h"

// constructor
NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::NodeCoordinatesDisplacementsVelocitiesAccelerations2D()
:NodeCoordinates2D(),NodeDisplacements2D(),NodeVelocities2D(),NodeAccelerations2D()
{
}

// sets the global dofs
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::SetGlobalDofs(int& rDOF)
{
    NodeCoordinates2D::SetGlobalDofs(rDOF);
    NodeDisplacements2D::SetGlobalDofs(rDOF);
}

// write dof values to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    NodeDisplacements2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
}

// extract dof values from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    NodeDisplacements2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
}

// renumber the global dofs according to predefined ordering
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    NodeCoordinates2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    NodeDisplacements2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
}

// write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // set velocities based on the displacement dofs
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 2; dim_count++)
    {
		int dof = this->NodeDisplacements2D::mDOF[dim_count];
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
		this->NodeVelocities2D::mVelocities[dim_count] = value;
    }
}

// extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // based on displacement dof number write velocities into global vectors
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 2; dim_count++)
    {
		int dof = this->NodeDisplacements2D::mDOF[dim_count];
		assert(dof >= 0);
		double value = this->NodeVelocities2D::mVelocities[dim_count];
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
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    NodeCoordinates2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // set accelerations based on the displacement dofs
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 2; dim_count++)
    {
		int dof = this->NodeDisplacements2D::mDOF[dim_count];
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
		this->NodeAccelerations2D::mAccelerations[dim_count] = value;
    }
}

// extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
void NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    NodeCoordinates2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);

    // based on displacement dof number write accelerations into global vectors
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for(int dim_count=0; dim_count < 2; dim_count++)
    {
		int dof = this->NodeDisplacements2D::mDOF[dim_count];
		assert(dof >= 0);
		double value = this->NodeAccelerations2D::mAccelerations[dim_count];
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
std::string NuTo::NodeCoordinatesDisplacementsVelocitiesAccelerations2D::GetNodeTypeStr()const
{
	return std::string("NodeCoordinatesDisplacementsVelocitiesAccelerations2D");
}
