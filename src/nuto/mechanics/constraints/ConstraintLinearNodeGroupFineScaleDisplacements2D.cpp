// $Id: ConstraintLinearNodeGroupFineScaleDisplacements2D.cpp 314 2010-09-27 16:31:43Z unger3 $
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeDisplacementsMultiscale2D.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupFineScaleDisplacements2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::ConstraintLinearNodeGroupFineScaleDisplacements2D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue) :
        ConstraintNodeGroup(rGroup), ConstraintLinear()
{
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::ConstraintLinearNodeGroupFineScaleDisplacements2D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.mEigenMatrix.data(),2*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::ConstraintLinearNodeGroupFineScaleDisplacements2D] direction vector has zero length.");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mRHS = rValue;

    mConstrainToStructureMultiscalePeriodicDofs = false;
    mStructure = 0;
}

//! @brief constructor
NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::ConstraintLinearNodeGroupFineScaleDisplacements2D(const StructureMultiscale* rStructure, const Group<NodeBase>* rGroup) :
        ConstraintNodeGroup(rGroup), ConstraintLinear()
{
	mStructure =  rStructure;
	mRHS = 0;
    mDirection[0]=0;
    mDirection[1]=0;
    mConstrainToStructureMultiscalePeriodicDofs = true;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::GetNumLinearConstraints()const
{
	if (mConstrainToStructureMultiscalePeriodicDofs)
    	return 2*mGroup->GetNumMembers();
    else
    	return mGroup->GetNumMembers();
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    if (mConstrainToStructureMultiscalePeriodicDofs)
    {
		for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
		{
			rRHS(curConstraintEquation,0) = 0;
			if (itNode->second->GetNumFineScaleDisplacements()!=2)
			{
				throw MechanicsException("[NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::AddToConstraintMatrix] Node does not have fine scale displacements or has more than two displacement components.");
			}
    		NodeDisplacementsMultiscale2D* nodePtr=itNode->second->AsNodeDisplacementsMultiscale2D();
    		rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[0],nodePtr->GetShapeFunctionMultiscalePeriodicX()[0]);
    		rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[1],nodePtr->GetShapeFunctionMultiscalePeriodicX()[1]);
    		rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[2],nodePtr->GetShapeFunctionMultiscalePeriodicX()[2]);
			rConstraintMatrix.AddEntry(curConstraintEquation,itNode->second->GetDofFineScaleDisplacement(0),-1);
    		//std::cout << "add constraint with " << itNode->second->GetDofFineScaleDisplacement(0) << " and " << mStructure->GetDOFPeriodicBoundaryDisplacements()[0] << " " << mStructure->GetDOFPeriodicBoundaryDisplacements()[1] << " "<<mStructure->GetDOFPeriodicBoundaryDisplacements()[2]<< std::endl;

    		curConstraintEquation++;

    		rRHS(curConstraintEquation,0) = 0;
    		rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[0],nodePtr->GetShapeFunctionMultiscalePeriodicY()[0]);
    		rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[1],nodePtr->GetShapeFunctionMultiscalePeriodicY()[1]);
    		rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[2],nodePtr->GetShapeFunctionMultiscalePeriodicY()[2]);
    		rConstraintMatrix.AddEntry(curConstraintEquation,itNode->second->GetDofFineScaleDisplacement(1),-1);
    		//std::cout << "add constraint with " << itNode->second->GetDofFineScaleDisplacement(1) << " and " << mStructure->GetDOFPeriodicBoundaryDisplacements()[0] << " " << mStructure->GetDOFPeriodicBoundaryDisplacements()[1] << " "<<mStructure->GetDOFPeriodicBoundaryDisplacements()[2]<< std::endl;

    		curConstraintEquation++;
		}
    }
    else
    {
    	for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
		{
			rRHS(curConstraintEquation,0) = mRHS;
			if (itNode->second->GetNumFineScaleDisplacements()!=2)
			{
				throw MechanicsException("[NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::AddToConstraintMatrix] Node does not have fine scale displacements or has more than two displacement components.");
			}
			if (fabs(mDirection[0])>1e-18)
			{
				rConstraintMatrix.AddEntry(curConstraintEquation,itNode->second->GetDofFineScaleDisplacement(0),mDirection[0]);
			}
			if (fabs(mDirection[1])>1e-18)
			{
				rConstraintMatrix.AddEntry(curConstraintEquation,itNode->second->GetDofFineScaleDisplacement(1),mDirection[1]);
			}

			curConstraintEquation++;
		}
    }
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeGroupFineScaleDisplacements2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNodeGroup)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mRHS)
       & BOOST_SERIALIZATION_NVP(mDirection)
       & BOOST_SERIALIZATION_NVP(mConstrainToStructureMultiscalePeriodicDofs);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeGroupFineScaleDisplacements2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D)
#endif // ENABLE_SERIALIZATION
