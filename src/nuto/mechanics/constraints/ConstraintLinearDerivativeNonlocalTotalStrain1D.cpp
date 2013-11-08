// $Id: ConstraintLinearDerivativeNonlocalTotalStrain1D.cpp 625 2013-04-22 16:37:11Z unger3 $

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
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintLinearDerivativeNonlocalTotalStrain1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/elements/Truss.h"
#include "nuto/mechanics/elements/BoundaryGradientDamage1D.h"

// constructor
NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::ConstraintLinearDerivativeNonlocalTotalStrain1D(const ElementBase* rParentElement, double rLocalIpCoordinate):
        ConstraintLinear()
{
    // set parent element
	mParentElement = rParentElement;
	mLocalIpCoordinate = rLocalIpCoordinate;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::SetRHS(double rRHS)
{
	throw MechanicsException("[NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::SetRHS] the right hand side is always zero.");
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix)const
{
	switch (mParentElement->GetEnumType())
	{
	case Element::TRUSS1D2N:
	case Element::TRUSS1D3N:
	{
		const Truss* elementPtr (mParentElement->AsTruss());
		// add constraint to constrain matrix
		int numNonlocalTotalStrain(elementPtr->GetNumShapeFunctionsNonlocalTotalStrain());

		std::vector<double> derivativeShapeFunctionsNaturalNonlocalTotalStrain(numNonlocalTotalStrain);  //allocate space for derivatives of shape functions
		//std::vector<double> derivativeShapeFunctionsLocalNonlocalTotalStrain(numNonlocalTotalStrain);    //allocate space for derivatives of shape functions

		//derivative in natural coordinate system
		elementPtr->CalculateDerivativeShapeFunctionsNonlocalTotalStrain(mLocalIpCoordinate, derivativeShapeFunctionsNaturalNonlocalTotalStrain);

		//derivative in local coordinate system
		//for (unsigned int count=0; count<derivativeShapeFunctionsLocalNonlocalTotalStrain.size(); count++)
		//{
		//	derivativeShapeFunctionsLocalNonlocalTotalStrain[count] = derivativeShapeFunctionsNaturalNonlocalTotalStrain[count]/detJ;
		//}
		// For 1D, there is only one point, so the detJ can be neglected

		for (int count=0; count<numNonlocalTotalStrain; count++)
		{
			rConstraintMatrix.AddValue(curConstraintEquation,
					elementPtr->GetNodeNonlocalTotalStrain(count)->GetDofNonlocalTotalStrain(0),
	        		derivativeShapeFunctionsNaturalNonlocalTotalStrain[count]);
		}
	}
	break;
	case Element::BOUNDARYGRADIENTDAMAGE1D:
	{
		const BoundaryGradientDamage1D* elementPtr (mParentElement->AsBoundaryGradientDamage1D());
		// add constraint to constrain matrix
		int numNonlocalTotalStrain(elementPtr->GetNumNodesField());

		std::vector<double> derivativeShapeFunctionsNaturalNonlocalTotalStrain(numNonlocalTotalStrain);  //allocate space for derivatives of shape functions
		//std::vector<double> derivativeShapeFunctionsLocalNonlocalTotalStrain(numNonlocalTotalStrain);    //allocate space for derivatives of shape functions

		//derivative in natural coordinate system
		elementPtr->CalculateDerivativeShapeFunctionsField(mLocalIpCoordinate, derivativeShapeFunctionsNaturalNonlocalTotalStrain);

		//derivative in local coordinate system
		//for (unsigned int count=0; count<derivativeShapeFunctionsLocalNonlocalTotalStrain.size(); count++)
		//{
		//	derivativeShapeFunctionsLocalNonlocalTotalStrain[count] = derivativeShapeFunctionsNaturalNonlocalTotalStrain[count]/detJ;
		//}
		// For 1D, there is only one point, so the detJ can be neglected

		for (int count=0; count<numNonlocalTotalStrain; count++)
		{
			rConstraintMatrix.AddValue(curConstraintEquation,
					elementPtr->GetNode(count)->GetDofNonlocalTotalStrain(0),
	        		derivativeShapeFunctionsNaturalNonlocalTotalStrain[count]);
		}
	}
	break;
	default:
		throw MechanicsException("[NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::AddToConstraintMatrix] unsopported element type.");
	}



    // increase constraint equation number
    curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const
{
    // set right hand side value
    rRHS(curConstraintEquation) = 0;

    // increase constraint equation number
    curConstraintEquation++;
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearDerivativeNonlocalTotalStrain1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mParentElement)
       & BOOST_SERIALIZATION_NVP(mLocalIpCoordinate);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearDerivativeNonlocalTotalStrain1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D)
#endif // ENABLE_SERIALIZATION
