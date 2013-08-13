// $Id$

#ifndef NODE_DOF_H
#define NODE_DOF_H
#include "nuto/mechanics/nodes/NodeDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures, TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain >::NodeDof() //: NuTo::NodeBase()
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
		for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
		{
			mNonlocalEqPlasticStrain[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumNonlocalTotalStrain; count++)
		{
			mNonlocalTotalStrain[countDerivatives][count]=0.;
		}
	}
	for (int count=0; count<TNumDisplacements; count++)
	{
		mDofDisplacements[count]=-1;
	}
	for (int count=0; count<TNumRotations; count++)
	{
		mDofRotations[count]=-1;
	}
	for (int count=0; count<TNumTemperatures; count++)
	{
		mDofTemperatures[count]=-1;
	}
	for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
	{
		mDofNonlocalEqPlasticStrain[count]=-1;
	}
	for (int count=0; count<TNumNonlocalTotalStrain; count++)
	{
		mDofNonlocalTotalStrain[count]=-1;
	}

}

//! @brief ... destructor
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::~NodeDof() //: NuTo::~NodeBase()
{
}

//! @brief assignment operator
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
     operator=(NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain> const& rOther)
{
    mDisplacements = rOther.mDisplacements;
    mDofDisplacements = rOther.mDofDisplacements;

    mRotations = rOther.mRotations;
    mDofRotations = rOther.mDofRotations;

    mTemperatures = rOther.mTemperatures;
    mDofTemperatures = rOther.mDofTemperatures;

    mNonlocalEqPlasticStrain = rOther.mNonlocalEqPlasticStrain;
    mDofNonlocalEqPlasticStrain = rOther.mDofNonlocalEqPlasticStrain;

    mNonlocalTotalStrain = rOther.mNonlocalTotalStrain;
    mDofNonlocalTotalStrain = rOther.mDofNonlocalTotalStrain;
}

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
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
    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
    	mDofNonlocalEqPlasticStrain[count]=rDOF++;
    }
    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
    	mDofNonlocalTotalStrain[count]=rDOF++;
    }
}

//! @brief write dof values to the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures, TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain>::
SetGlobalDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues)
{
    assert(TNumTimeDerivatives>=rTimeDerivative);

    for (int count=0; count<TNumDisplacements; count++)
    {
        int dof = this->mDofDisplacements[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
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
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
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
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mTemperatures[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
        int dof = this->mDofNonlocalEqPlasticStrain[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        if (value<0)
        {
        	std::cout << "value is less then zero" << value << std::endl;
        }
        this->mNonlocalEqPlasticStrain[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
        int dof = this->mDofNonlocalTotalStrain[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mNonlocalTotalStrain[rTimeDerivative][count] = value;
    }
}


//! @brief extract dof values from the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetGlobalDofValues(int rTimeDerivative, NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
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
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
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
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
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
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }

    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
        int dof = this->mDofNonlocalEqPlasticStrain[count];
        double value = this->mNonlocalEqPlasticStrain[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }

    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
        int dof = this->mDofNonlocalTotalStrain[count];
        double value = this->mNonlocalTotalStrain[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }
}


//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
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

    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
    	mDofNonlocalEqPlasticStrain[count]=rMappingInitialToNewOrdering[mDofNonlocalEqPlasticStrain[count]];
    }

    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
    	mDofNonlocalTotalStrain[count]=rMappingInitialToNewOrdering[mDofNonlocalTotalStrain[count]];
    }
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNumTimeDerivatives()const
{
	return TNumTimeDerivatives;
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNumDisplacements()const
{
	return TNumDisplacements;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDofDisplacement(int rComponent)const
{
	assert(rComponent<TNumDisplacements);
	return mDofDisplacements[rComponent];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDisplacements1D(double rDisplacements[1])const
{
	assert(TNumDisplacements==1);
	rDisplacements[0] = mDisplacements[0][0];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDisplacements1D(int rTimeDerivative, double rDisplacements[1])const
{
	assert(TNumDisplacements==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetDisplacements1D(const double rDisplacements[1])
{
	assert(TNumDisplacements==1);
	mDisplacements[0][0] = rDisplacements[0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetDisplacements1D(int rTimeDerivative, const double rDisplacements[1])
{
	assert(TNumDisplacements==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDisplacements2D(double rDisplacements[2])const
{
	assert(TNumDisplacements==2);
	rDisplacements[0] = mDisplacements[0][0];
	rDisplacements[1] = mDisplacements[0][1];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDisplacements2D(int rTimeDerivative, double rDisplacements[2])const
{
	assert(TNumDisplacements==2);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
	rDisplacements[1] = mDisplacements[rTimeDerivative][1];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetDisplacements2D(const double rDisplacements[2])
{
	assert(TNumDisplacements==2);
	mDisplacements[0][0] = rDisplacements[0];
	mDisplacements[0][1] = rDisplacements[1];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetDisplacements2D(int rTimeDerivative, const double rDisplacements[2])
{
	assert(TNumDisplacements==2);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
	mDisplacements[rTimeDerivative][1] = rDisplacements[1];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDisplacements3D(double rDisplacements[3])const
{
	assert(TNumDisplacements==3);
	rDisplacements[0] = mDisplacements[0][0];
	rDisplacements[1] = mDisplacements[0][1];
	rDisplacements[2] = mDisplacements[0][2];
}

//! @brief returns the displacements of the node
//! @return displacement
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
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
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetDisplacements3D(const double rDisplacements[3])
{
	assert(TNumDisplacements==3);
	mDisplacements[0][0] = rDisplacements[0];
	mDisplacements[0][1] = rDisplacements[1];
	mDisplacements[0][2] = rDisplacements[2];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
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
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDisplacement(short rIndex)const
{
	assert(rIndex<TNumDisplacements);
	return mDisplacements[0][rIndex];
}


//! @brief returns the number of Rotations of the node
//! @return number of Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNumRotations()const
{
	return TNumRotations;
}

//! @brief gives the global DOF of a Rotation component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDofRotation(int rComponent)const
{
	assert(TNumRotations>rComponent);
	return mDofRotations[rComponent];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetRotations2D(double rRotations[1])const
{
	assert(TNumRotations==1);
	rRotations[0] = mRotations[0][0];
}

//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetRotations2D(int rTimeDerivative, double rRotations[1])const
{
	assert(TNumRotations==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rRotations[rTimeDerivative] = mRotations[0][0];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetRotations2D(const double rRotations[1])
{
	assert(TNumRotations==1);
	mRotations[0][0] = rRotations[0];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetRotations2D(int rTimeDerivative, const double rRotations[1])
{
	assert(TNumRotations==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mRotations[rTimeDerivative][0] = rRotations[0];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetRotations3D(double rRotations[3])const
{
	assert(TNumRotations==3);
	rRotations[0] = mRotations[0][0];
	rRotations[1] = mRotations[0][1];
	rRotations[2] = mRotations[0][2];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetRotations3D(int rTimeDerivative, double rRotations[3])const
{
	assert(TNumRotations==3);
	rRotations[0] = mRotations[rTimeDerivative][0];
	rRotations[1] = mRotations[rTimeDerivative][1];
	rRotations[2] = mRotations[rTimeDerivative][2];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetRotations3D(const double rRotations[3])
{
	assert(TNumRotations==3);
	mRotations[0][0] = rRotations[0];
	mRotations[0][1] = rRotations[1];
	mRotations[0][2] = rRotations[2];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetRotations3D(int rTimeDerivative, const double rRotations[3])
{
	assert(TNumRotations==3);
	mRotations[rTimeDerivative][0] = rRotations[0];
	mRotations[rTimeDerivative][1] = rRotations[1];
	mRotations[rTimeDerivative][2] = rRotations[2];
}


//! @brief returns the Rotations of the node
//! @return Rotation
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetRotation(short rIndex)const
{
	assert(TNumRotations>rIndex);
	return mRotations[0][rIndex];
}

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNumTemperatures()const
{
	return TNumTemperatures;
}

//! @brief returns the temperature of the node
//! @return temperature
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetTemperature()const
{
	assert(TNumTemperatures==1);
	return mTemperatures[0][0];
}

//! @brief returns the temperature of the node
//! @return temperature
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetTemperature(int rTimeDerivative)const
{
	assert(TNumTemperatures==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	return mTemperatures[rTimeDerivative][0];
}

//! @brief set the temperature of the node
//! @param rTemperature  given temperature
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetTemperature(double rTemperatures)
{
	assert(TNumTemperatures==1);
	mTemperatures[0][0] = rTemperatures;
}

//! @brief set the temperature of the node
//! @param rTemperature  given temperature
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetTemperature(int rTimeDerivative, double rTemperatures)
{
	assert(TNumTemperatures==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mTemperatures[rTimeDerivative][0] = rTemperatures;
}

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDofTemperature()const
{
	assert(TNumTemperatures==1);
	return mDofTemperatures[0];

}

//! @brief returns the number of Damage of the node
//! @return number of Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNumNonlocalEqPlasticStrain()const
{
	return TNumNonlocalEqPlasticStrain;
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalEqPlasticStrain(double* rNonlocalEqPlasticStrain)const
{
	assert(TNumNonlocalEqPlasticStrain==2);
	rNonlocalEqPlasticStrain[0] = mNonlocalEqPlasticStrain[0][0];
	rNonlocalEqPlasticStrain[1] = mNonlocalEqPlasticStrain[0][1];
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalEqPlasticStrain(int rTimeDerivative, double* rNonlocalEqPlasticStrain)const
{
	assert(TNumNonlocalEqPlasticStrain==2);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalEqPlasticStrain[0] = mNonlocalEqPlasticStrain[rTimeDerivative][0];
	rNonlocalEqPlasticStrain[1] = mNonlocalEqPlasticStrain[rTimeDerivative][1];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalEqPlasticStrain(const double* rNonlocalEqPlasticStrain)
{
	assert(TNumNonlocalEqPlasticStrain==2);
	mNonlocalEqPlasticStrain[0][0] = rNonlocalEqPlasticStrain[0];
	mNonlocalEqPlasticStrain[0][1] = rNonlocalEqPlasticStrain[1];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumNonlocalEqPlasticStrains, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumNonlocalEqPlasticStrains,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalEqPlasticStrain(int rTimeDerivative, const double* rNonlocalEqPlasticStrain)
{
	assert(TNumNonlocalEqPlasticStrain==2);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalEqPlasticStrain[rTimeDerivative][0] = rNonlocalEqPlasticStrain[0];
	mNonlocalEqPlasticStrain[rTimeDerivative][1] = rNonlocalEqPlasticStrain[1];
}

//! @brief gives the global DOF of a Damage component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDofNonlocalEqPlasticStrain(int rComponent)const
{
	assert(rComponent<TNumNonlocalEqPlasticStrain);
	return mDofNonlocalEqPlasticStrain[rComponent];

}

//! @brief returns the number of Damage of the node
//! @return number of Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNumNonlocalTotalStrain()const
{
	return TNumNonlocalTotalStrain;
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain1D(double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==1);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[0][0];
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain2D(double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==3);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[0][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[0][1];
	rNonlocalTotalStrain[2] = mNonlocalTotalStrain[0][2];
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain3D(double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==6);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[0][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[0][1];
	rNonlocalTotalStrain[2] = mNonlocalTotalStrain[0][2];
	rNonlocalTotalStrain[3] = mNonlocalTotalStrain[0][3];
	rNonlocalTotalStrain[4] = mNonlocalTotalStrain[0][4];
	rNonlocalTotalStrain[5] = mNonlocalTotalStrain[0][5];
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain1D(int rTimeDerivative, double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[rTimeDerivative][0];
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain2D(int rTimeDerivative, double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==3);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[rTimeDerivative][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[rTimeDerivative][1];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[rTimeDerivative][2];
}

//! @brief returns the Damage of the node
//! @return Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain3D(int rTimeDerivative, double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==6);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[rTimeDerivative][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[rTimeDerivative][1];
	rNonlocalTotalStrain[2] = mNonlocalTotalStrain[rTimeDerivative][2];
	rNonlocalTotalStrain[3] = mNonlocalTotalStrain[rTimeDerivative][3];
	rNonlocalTotalStrain[4] = mNonlocalTotalStrain[rTimeDerivative][4];
	rNonlocalTotalStrain[5] = mNonlocalTotalStrain[rTimeDerivative][5];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalTotalStrain1D(const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==1);
	mNonlocalTotalStrain[0][0] = rNonlocalTotalStrain[0];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalTotalStrain2D(const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==3);
	mNonlocalTotalStrain[0][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[0][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[0][2] = rNonlocalTotalStrain[2];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalTotalStrain3D(const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==6);
	mNonlocalTotalStrain[0][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[0][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[0][2] = rNonlocalTotalStrain[2];
	mNonlocalTotalStrain[0][3] = rNonlocalTotalStrain[3];
	mNonlocalTotalStrain[0][4] = rNonlocalTotalStrain[4];
	mNonlocalTotalStrain[0][5] = rNonlocalTotalStrain[5];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalTotalStrain1D(int rTimeDerivative, const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalTotalStrain[rTimeDerivative][0] = rNonlocalTotalStrain[0];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalTotalStrain2D(int rTimeDerivative, const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==3);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalTotalStrain[rTimeDerivative][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[rTimeDerivative][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[rTimeDerivative][2] = rNonlocalTotalStrain[2];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
SetNonlocalTotalStrain3D(int rTimeDerivative, const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==6);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalTotalStrain[rTimeDerivative][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[rTimeDerivative][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[rTimeDerivative][2] = rNonlocalTotalStrain[2];
	mNonlocalTotalStrain[rTimeDerivative][3] = rNonlocalTotalStrain[3];
	mNonlocalTotalStrain[rTimeDerivative][4] = rNonlocalTotalStrain[4];
	mNonlocalTotalStrain[rTimeDerivative][5] = rNonlocalTotalStrain[5];
}

//! @brief returns the nonlocal total strain component of the node
//! @return strain component (rTimeDerivative=0)
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
double NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetNonlocalTotalStrain(short rIndex)const
{
	assert(rIndex<TNumNonlocalTotalStrain);
	return mNonlocalTotalStrain[0][rIndex];
}

//! @brief gives the global DOF of a Damage component
//! @param rComponent component
//! @return global DOF
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
int NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
GetDofNonlocalTotalStrain(int rComponent)const
{
	assert(rComponent<TNumNonlocalTotalStrain);
	return mDofNonlocalTotalStrain[rComponent];

}

//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
std::string NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
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
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
void NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{

}
#endif // ENABLE_VISUALIZE

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain >
NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain >* NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>::
Clone()const
{
    return new NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures,TNumNonlocalEqPlasticStrain,TNumNonlocalTotalStrain>(*this);
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

