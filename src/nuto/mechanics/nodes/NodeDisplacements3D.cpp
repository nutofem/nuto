#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeDisplacements3D.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
NuTo::NodeDisplacements3D::NodeDisplacements3D() : NodeBase ()
{
	this->mDisplacements[0]=0;
	this->mDisplacements[1]=0;
	this->mDisplacements[2]=0;
    this->mDOF[0]=-1;
    this->mDOF[1]=-1;
    this->mDOF[2]=-1;
}

//! @brief constructor
NuTo::NodeDisplacements3D::NodeDisplacements3D (const double rDisplacements[3])  : NodeBase ()
{
	this->mDisplacements[0]= rDisplacements[0];
	this->mDisplacements[1]= rDisplacements[1];
	this->mDisplacements[2]= rDisplacements[2];
    this->mDOF[0]=-1;
    this->mDOF[1]=-1;
    this->mDOF[2]=-1;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeDisplacements3D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mDisplacements);
}
#endif // ENABLE_SERIALIZATION

//! @brief returns the number of displacements of the node
//! @return number of displacements
int NuTo::NodeDisplacements3D::GetNumDisplacements()const
{
    return 3;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeDisplacements3D::GetDofDisplacement(int rComponent)const
{
	if (rComponent<0 || rComponent>2)
		throw MechanicsException("[NuTo::NodeDisplacements2D::GetDofDisplacement] Node has only two displacement components.");
    return mDOF[rComponent];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeDisplacements3D::SetDisplacements3D(const double rDisplacements[1])
{
    mDisplacements[0] = rDisplacements[0];
    mDisplacements[1] = rDisplacements[1];
    mDisplacements[2] = rDisplacements[2];
}

//! @brief writes the displacements of a node to the prescribed pointer
//! @param rDisplacements displacements
void NuTo::NodeDisplacements3D::GetDisplacements3D(double rDisplacements[1])const
{
    rDisplacements[0] = mDisplacements[0];
    rDisplacements[1] = mDisplacements[1];
    rDisplacements[2] = mDisplacements[2];
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeDisplacements3D::SetGlobalDofs(int& rDOF)
{
    mDOF[0]=rDOF++;
    mDOF[1]=rDOF++;
    mDOF[2]=rDOF++;
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements3D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for (int count=0; count<3; count++)
    {
        int dof = this->mDOF[count];
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
        this->mDisplacements[count] = value;
    }
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements3D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for (int count=0; count<3; count++)
    {
        int dof = this->mDOF[count];
        double value = this->mDisplacements[count];
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

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements3D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeDisplacements3D::SetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements3D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeDisplacements3D::GetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements3D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeDisplacements3D::SetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements3D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeDisplacements3D::GetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeDisplacements3D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    mDOF[0]=rMappingInitialToNewOrdering[mDOF[0]];
    mDOF[1]=rMappingInitialToNewOrdering[mDOF[1]];
    mDOF[2]=rMappingInitialToNewOrdering[mDOF[2]];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeDisplacements3D::GetNodeTypeStr()const
{
	return std::string("NodeDisplacements3D");
}

