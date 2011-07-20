// $Id: NodeDisplacements2D.cpp 343 2010-10-19 07:43:10Z arnold2 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeDisplacementsMultiscale2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"

//! @brief constructor
NuTo::NodeDisplacementsMultiscale2D::NodeDisplacementsMultiscale2D(NuTo::StructureMultiscale* rStructureMultiscale, bool rCrackedDomain) : NodeBase ()
{
	this->mCrackedDomain = rCrackedDomain;
	this->mFineScaleDisplacements[0]=0.;
    this->mFineScaleDisplacements[1]=0.;
    this->mDOF[0]=-1;
    this->mDOF[1]=-1;
    mStructureMultiscale = rStructureMultiscale;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeDisplacementsMultiscale2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacementsMultiscale2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacementsMultiscale2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacementsMultiscale2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacementsMultiscale2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacementsMultiscale2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeDisplacementsMultiscale2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeDisplacementsMultiscale2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mFineScaleDisplacements)
       & BOOST_SERIALIZATION_NVP(mShapeFunctionX)
       & BOOST_SERIALIZATION_NVP(mShapeFunctionY)
       & BOOST_SERIALIZATION_NVP(mCrackedDomain)
       & BOOST_SERIALIZATION_NVP(mDOF)
       & BOOST_SERIALIZATION_NVP(mStructureMultiscale);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeDisplacementsMultiscale2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeDisplacementsMultiscale2D)
BOOST_CLASS_TRACKING(NuTo::NodeDisplacementsMultiscale2D, track_always)
#endif // ENABLE_SERIALIZATION

//! @brief returns the number of displacements of the node
//! @return number of displacements
int NuTo::NodeDisplacementsMultiscale2D::GetNumFineScaleDisplacements()const
{
    return 2;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeDisplacementsMultiscale2D::GetDofFineScaleDisplacement(int rComponent)const
{
    if (rComponent<0 || rComponent>1)
        throw MechanicsException("[NuTo::NodeDisplacementsMultiscale2D::GetDofFineScaleDisplacement] Node has only two displacement components.");
    return mDOF[rComponent];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeDisplacementsMultiscale2D::SetFineScaleDisplacements2D(const double rDisplacements[2])
{
    mFineScaleDisplacements[0] = rDisplacements[0];
    mFineScaleDisplacements[1] = rDisplacements[1];
}

//! @brief writes the displacements of a node to the prescribed pointer
//! @param rDisplacements displacements
void NuTo::NodeDisplacementsMultiscale2D::GetFineScaleDisplacements2D(double rDisplacements[2])const
{
    rDisplacements[0] = mFineScaleDisplacements[0];
    rDisplacements[1] = mFineScaleDisplacements[1];
}

//! @brief writes the displacements of a node to the prescribed pointer
//! the difference is e.g. using XFEM, when the nodal degrees of freedom are not identical
//! @param rDisplacements displacements
void NuTo::NodeDisplacementsMultiscale2D::GetDisplacements2D(double rDisplacements[2])const
{
    double coordinates[2], displacements[2];
    GetCoordinates2D(coordinates);
    if (mCrackedDomain)
    {
    	mStructureMultiscale->GetDisplacementsEpsilonHom2D(coordinates, rDisplacements, mStructureMultiscale->GetCenterDamage());
    	mStructureMultiscale->GetDisplacementsCrack2D(coordinates, displacements);
	    rDisplacements[0] += displacements[0] + mFineScaleDisplacements[0];
	    rDisplacements[1] += displacements[1] + mFineScaleDisplacements[1];
    }
    else
    {
		mStructureMultiscale->GetDisplacementsEpsilonHom2D(coordinates, rDisplacements, mStructureMultiscale->GetCenterHomogeneous());
	    rDisplacements[0] += mFineScaleDisplacements[0];
	    rDisplacements[1] += mFineScaleDisplacements[1];
    }
}

//! @brief writes the displacements of a node to the prescribed pointer
//! the difference is e.g. using XFEM, when the nodal degrees of freedom are not identical
//! @param rDisplacements displacements
double NuTo::NodeDisplacementsMultiscale2D::GetDisplacement(short rIndex)const
{
    double coordinates[2], displacements[2],displacementsCrack[2];
    GetCoordinates2D(coordinates);

    if (mCrackedDomain)
    {
		mStructureMultiscale->GetDisplacementsEpsilonHom2D(coordinates, displacements, mStructureMultiscale->GetCenterDamage());
		mStructureMultiscale->GetDisplacementsCrack2D(coordinates, displacementsCrack);
	    displacements[0] += displacementsCrack[0] + mFineScaleDisplacements[0];
	    displacements[1] += displacementsCrack[1] + mFineScaleDisplacements[1];
    }
    else
    {
		mStructureMultiscale->GetDisplacementsEpsilonHom2D(coordinates, displacements, mStructureMultiscale->GetCenterHomogeneous());
	    displacements[0] += mFineScaleDisplacements[0];
	    displacements[1] += mFineScaleDisplacements[1];
    }

    if (rIndex==0 || rIndex==1)
        return displacements[rIndex];
    else
        throw MechanicsException("[NuTo::NodeDisplacementsMultiscale2D::GetDisplacement] node has only two displacements.");
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeDisplacementsMultiscale2D::SetGlobalDofs(int& rDOF)
{
    mDOF[0]=rDOF++;
    mDOF[1]=rDOF++;
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacementsMultiscale2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for (int count=0; count<2; count++)
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
        this->mFineScaleDisplacements[count] = value;
    }
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacementsMultiscale2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    for (int count=0; count<2; count++)
    {
        int dof = this->mDOF[count];
        double value = this->mFineScaleDisplacements[count];
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
void NuTo::NodeDisplacementsMultiscale2D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    throw MechanicsException("[NuTo::NodeDisplacementsMultiscale2D::SetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacementsMultiscale2D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    throw MechanicsException("[NuTo::NodeDisplacementsMultiscale2D::GetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacementsMultiscale2D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    throw MechanicsException("[NuTo::NodeDisplacementsMultiscale2D::SetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacementsMultiscale2D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    throw MechanicsException("[NuTo::NodeDisplacementsMultiscale2D::GetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeDisplacementsMultiscale2D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    mDOF[0]=rMappingInitialToNewOrdering[mDOF[0]];
    mDOF[1]=rMappingInitialToNewOrdering[mDOF[1]];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeDisplacementsMultiscale2D::GetNodeTypeStr()const
{
    return std::string("NodeDisplacementsMultiscale2D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeDisplacementsMultiscale2D::GetNodeType()const
{
    return Node::NodeDisplacementsMultiscale2D;
}

//! @brief scales the shape functions
//! @parameter rShapeFunction  (1..3 corresponding to macro strains exx, eyy, and gxy)
//! @parameter rScalingFactor rScalingFactor
void NuTo::NodeDisplacementsMultiscale2D::ScaleShapeFunctionMultiscalePeriodic(int rShapeFunction, double rScalingFactor)
{
	assert(rShapeFunction>=0 && rShapeFunction<=2);
	mShapeFunctionX[rShapeFunction]*=rScalingFactor;
	mShapeFunctionY[rShapeFunction]*=rScalingFactor;
}

//! @brief set the shape functions based on the actual oscillations
//! @parameter shape function number (0..2)
void NuTo::NodeDisplacementsMultiscale2D::SetShapeFunctionMultiscalePeriodic(int rShapeFunction)
{
	assert(rShapeFunction>-1 && rShapeFunction<3);
	mShapeFunctionX[rShapeFunction] = mFineScaleDisplacements[0];
	mShapeFunctionY[rShapeFunction] = mFineScaleDisplacements[1];
}
