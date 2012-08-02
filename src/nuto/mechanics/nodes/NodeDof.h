// $Id$

#ifndef NODE_DOF_H
#define NODE_DOF_H
#include "nuto/mechanics/nodes/NodeDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::NodeDof() //: NuTo::NodeBase()
{
	for (int countDerivatives=0; countDerivatives<TNumTimeDerivatives+1; countDerivatives++)
	{
		for (int count=0; count<TNumDisplacements; count++)
		{
			mDisplacements[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumRotations; count++)
		{
			mRotations[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumTemperatures; count++)
		{
			mTemperatures[countDerivatives][count]=0.;
		}
	}
}

//! @brief ... destructor
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::~NodeDof() //: NuTo::~NodeBase()
{
}

//! @brief assignment operator
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
     operator=(NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures> const& rOther)
{
    mDisplacements = rOther.mDisplacements;
    mDofDisplacements = rOther.mDofDisplacements;

    mRotations = rOther.mRotations;
    mDofRotations = rOther.mDofRotations;

    mTemperatures = rOther.mTemperatures;
    mDofTemperatures = rOther.mDofTemperatures;
}

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetGlobalDofs(int& rDOF)
{
    for (int count=0; count<TNumDisplacements; count++)
    {
    	mDofDisplacements[count]=rDOF++;
    }
    for (int count=0; count<TNumRotations; count++)
    {
    	mDofRotations[count]=rDOF++;
    }
    for (int count=0; count<TNumTemperatures; count++)
    {
    	mDofTemperatures[count]=rDOF++;
    }
}

//! @brief write dof values to the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetGlobalDofValues(int rTimeDerivative, const NuTo::FullMatrix<double>& rActiveDofValues, const NuTo::FullMatrix<double>& rDependentDofValues)
{
	assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    assert(TNumTimeDerivatives>=rTimeDerivative);

    for (int count=0; count<TNumDisplacements; count++)
    {
        int dof = this->mDofDisplacements[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof,0);
        }
        else
        {
            value = rActiveDofValues(dof,0);
        }
        this->mDisplacements[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumRotations; count++)
    {
        int dof = this->mDofRotations[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof,0);
        }
        else
        {
            value = rActiveDofValues(dof,0);
        }
        this->mRotations[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumTemperatures; count++)
    {
        int dof = this->mDofTemperatures[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof,0);
        }
        else
        {
            value = rActiveDofValues(dof,0);
        }
        this->mTemperatures[rTimeDerivative][count] = value;
    }
}


//! @brief extract dof values from the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetGlobalDofValues(int rTimeDerivative, NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const
{
	assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    assert(TNumTimeDerivatives>=rTimeDerivative);

    for (int count=0; count<TNumDisplacements; count++)
    {
        int dof = this->mDofDisplacements[count];
        double value = this->mDisplacements[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof,0) = value;
        }
        else
        {
            rActiveDofValues(dof,0) = value;
        }
    }

    for (int count=0; count<TNumRotations; count++)
    {
        int dof = this->mDofRotations[count];
        double value = this->mRotations[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof,0) = value;
        }
        else
        {
            rActiveDofValues(dof,0) = value;
        }
    }

    for (int count=0; count<TNumTemperatures; count++)
    {
        int dof = this->mDofTemperatures[count];
        double value = this->mTemperatures[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof,0) = value;
        }
        else
        {
            rActiveDofValues(dof,0) = value;
        }
    }
}


//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    for (int count=0; count<TNumDisplacements; count++)
    {
    	mDofDisplacements[count]=rMappingInitialToNewOrdering[mDofDisplacements[count]];
    }

    for (int count=0; count<TNumRotations; count++)
    {
    	mDofRotations[count]=rMappingInitialToNewOrdering[mDofRotations[count]];
    }

    for (int count=0; count<TNumTemperatures; count++)
    {
    	mDofTemperatures[count]=rMappingInitialToNewOrdering[mDofTemperatures[count]];
    }
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNumTimeDerivatives()const
{
	return TNumTimeDerivatives;
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNumDisplacements()const
{
	return TNumDisplacements;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDofDisplacement(int rComponent)const
{
	assert(rComponent<TNumDisplacements);
	return mDofDisplacements[rComponent];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacements1D(double rDisplacements[1])const
{
	assert(TNumDisplacements==1);
	rDisplacements[0] = mDisplacements[0][0];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacements1D(int rTimeDerivative, double rDisplacements[1])const
{
	assert(TNumDisplacements==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetDisplacements1D(const double rDisplacements[1])
{
	assert(TNumDisplacements==1);
	mDisplacements[0][0] = rDisplacements[0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetDisplacements1D(int rTimeDerivative, const double rDisplacements[1])
{
	assert(TNumDisplacements==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacements2D(double rDisplacements[2])const
{
	assert(TNumDisplacements==2);
	rDisplacements[0] = mDisplacements[0][0];
	rDisplacements[1] = mDisplacements[0][1];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacements2D(int rTimeDerivative, double rDisplacements[2])const
{
	assert(TNumDisplacements==2);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
	rDisplacements[1] = mDisplacements[rTimeDerivative][1];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetDisplacements2D(const double rDisplacements[2])
{
	assert(TNumDisplacements==2);
	mDisplacements[0][0] = rDisplacements[0];
	mDisplacements[0][1] = rDisplacements[1];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetDisplacements2D(int rTimeDerivative, const double rDisplacements[2])
{
	assert(TNumDisplacements==2);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
	mDisplacements[rTimeDerivative][1] = rDisplacements[1];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacements3D(double rDisplacements[3])const
{
	assert(TNumDisplacements==3);
	rDisplacements[0] = mDisplacements[0][0];
	rDisplacements[1] = mDisplacements[0][1];
	rDisplacements[2] = mDisplacements[0][2];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacements3D(int rTimeDerivative, double rDisplacements[3])const
{
	assert(TNumDisplacements==3);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
	rDisplacements[1] = mDisplacements[rTimeDerivative][1];
	rDisplacements[2] = mDisplacements[rTimeDerivative][2];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetDisplacements3D(const double rDisplacements[3])
{
	assert(TNumDisplacements==3);
	mDisplacements[0][0] = rDisplacements[0];
	mDisplacements[0][1] = rDisplacements[1];
	mDisplacements[0][2] = rDisplacements[2];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetDisplacements3D(int rTimeDerivative, const double rDisplacements[3])
{
	assert(TNumDisplacements==3);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
	mDisplacements[rTimeDerivative][1] = rDisplacements[1];
	mDisplacements[rTimeDerivative][2] = rDisplacements[2];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDisplacement(short rIndex)const
{
	assert(rIndex<TNumDisplacements);
	return mDisplacements[0][rIndex];
}


//! @brief returns the number of Rotations of the node
//! @return number of Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNumRotations()const
{
	return TNumRotations;
}

//! @brief gives the global DOF of a Rotation component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDofRotation(int rComponent)const
{
	assert(TNumRotations>rComponent);
	return mDofRotations[rComponent];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetRotations2D(double rRotations[1])const
{
	assert(TNumRotations==1);
	rRotations[0] = mRotations[0][0];
}

//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetRotations2D(int rTimeDerivative, double rRotations[1])const
{
	assert(TNumRotations==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rRotations[rTimeDerivative] = mRotations[0][0];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetRotations2D(const double rRotations[1])
{
	assert(TNumRotations==1);
	mRotations[0][0] = rRotations[0];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetRotations2D(int rTimeDerivative, const double rRotations[1])
{
	assert(TNumRotations==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mRotations[rTimeDerivative][0] = rRotations[0];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetRotations3D(double rRotations[3])const
{
	assert(TNumRotations==3);
	rRotations[0] = mRotations[0][0];
	rRotations[1] = mRotations[0][1];
	rRotations[2] = mRotations[0][2];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetRotations3D(int rTimeDerivative, double rRotations[3])const
{
	assert(TNumRotations==3);
	rRotations[0] = mRotations[rTimeDerivative][0];
	rRotations[1] = mRotations[rTimeDerivative][1];
	rRotations[2] = mRotations[rTimeDerivative][2];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetRotations3D(const double rRotations[3])
{
	assert(TNumRotations==3);
	mRotations[0][0] = rRotations[0];
	mRotations[0][1] = rRotations[1];
	mRotations[0][2] = rRotations[2];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetRotations3D(int rTimeDerivative, const double rRotations[3])
{
	assert(TNumRotations==3);
	mRotations[rTimeDerivative][0] = rRotations[0];
	mRotations[rTimeDerivative][1] = rRotations[1];
	mRotations[rTimeDerivative][2] = rRotations[2];
}


//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetRotation(short rIndex)const
{
	assert(TNumRotations>rIndex);
	return mRotations[0][rIndex];
}

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNumTemperatures()const
{
	return TNumTemperatures;
}

//! @brief returns the temperature of the node
//! @return temperature
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetTemperature(double* rTemperatures)const
{
	assert(TNumTemperatures==1);
	rTemperatures[0] = mTemperatures[0][0];
}

//! @brief set the temperature of the node
//! @param rTemperature  given temperature
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
SetTemperature(const double* rTemperatures)
{
	assert(TNumTemperatures==1);
	mTemperatures[0][0] = rTemperatures[0];
}

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetDofTemperature()const
{
	assert(TNumTemperatures==1);
	return mDofTemperatures[0];

}

//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
std::string NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNodeTypeStr()const
{
	throw MechanicsException("");
    return std::string("[NuTo::NodeDof::GetNodeTypeStr]to be done");
}

/*
//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
NuTo::Node::eNodeType NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNodeType()const
{
	return Nodes::
    NodeCoordinates1D,
    NodeCoordinates2D,
    NodeCoordinates3D,
    NodeCoordinatesTimeDerivative_0_Displacements1D,
    NodeCoordinatesTimeDerivative_0_Displacements2D,
    NodeCoordinatesTimeDerivative_0_Displacements3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsMultiscale1D,
    NodeCoordinatesTimeDerivative_0_DisplacementsMultiscale2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsMultiscale3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsNonlocalData2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsNonlocalData3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRadius2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRadius3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotations2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotations3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotationsRadius2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotationsRadius3D,
    NodeCoordinatesTimeDerivative_2_DisplacementsRotations2D,
    NodeCoordinatesTimeDerivative_2_DisplacementsRotationsRadius2D,
    NodeCoordinatesTimeDerivative_2_Displacements1D,
    NodeCoordinatesTimeDerivative_2_Displacements2D,
    NodeCoordinatesTimeDerivative_2_Displacements3D,
    NodeCoordinatesTimeDerivative_0_Temperature1D,
    NodeCoordinatesTimeDerivative_0_Temperature2D,
    NodeCoordinatesTimeDerivative_0_Temperature3D,

}
*/

#ifdef ENABLE_VISUALIZE
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{

}
#endif // ENABLE_VISUALIZE

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>* NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
Clone()const
{
    return new NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>(*this);
}


/*
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,1,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,1,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,1,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,3,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,1,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,3,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,0,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1,0,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,1,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,1,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,0,1,0>)))
#endif // ENABLE_SERIALIZATION
*/
#endif //NODEDOF_H

