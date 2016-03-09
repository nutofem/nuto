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
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constraints/ConstraintLinearDerivativeNonlocalTotalStrain1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"


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
	case Element::ELEMENT1D:
	{
        Eigen::VectorXd localIpCoordinate(1);
        localIpCoordinate(0) = mLocalIpCoordinate;

        const auto& interpolationTypeNonlocalTotalStrain = mParentElement->GetInterpolationType()->Get(Node::NONLOCALTOTALSTRAIN);
		//derivative in natural coordinate system
		auto derivativeShapeFunctionsNaturalNonlocalTotalStrain = interpolationTypeNonlocalTotalStrain.CalculateDerivativeShapeFunctionsNatural(localIpCoordinate);

		// For 1D, there is only one point, so the detJ can be neglected
		for (int count=0; count<interpolationTypeNonlocalTotalStrain.GetNumDofs(); count++)
		{
		    const NodeBase* node = mParentElement->GetNode(count, Node::NONLOCALTOTALSTRAIN);
			rConstraintMatrix.AddValue(curConstraintEquation,
			        node->GetDofNonlocalTotalStrain(0),
	        		derivativeShapeFunctionsNaturalNonlocalTotalStrain(count,0));
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
void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::save(Archive & ar, const unsigned int version) const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearDerivativeNonlocalTotalStrain1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear);

    std::uintptr_t mParentElementAddress = reinterpret_cast<std::uintptr_t>(mParentElement);
    ar & boost::serialization::make_nvp("mParentElement", mParentElementAddress);

    ar & BOOST_SERIALIZATION_NVP(mLocalIpCoordinate);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearDerivativeNonlocalTotalStrain1D" << std::endl;
#endif
}

template<class Archive>
void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearDerivativeNonlocalTotalStrain1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear);

    std::uintptr_t mParentElementAddress;
    ar & boost::serialization::make_nvp("mParentElement", mParentElementAddress);
    mParentElement = reinterpret_cast<const ElementBase*>(mParentElementAddress);

    ar & BOOST_SERIALIZATION_NVP(mLocalIpCoordinate);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearDerivativeNonlocalTotalStrain1D" << std::endl;
#endif
}

void NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D::SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast)
{
    std::map<std::uintptr_t, std::uintptr_t>::const_iterator it = mElementMapCast.find(reinterpret_cast<std::uintptr_t>(mParentElement));
    if (it!=mElementMapCast.end())
    {
        ElementBase** temp = const_cast<ElementBase**>(&mParentElement);
        *temp = reinterpret_cast<ElementBase*>(it->second);
    }
    else
        throw MechanicsException("[NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D] The ElementBase-Pointer could not be updated.");
}


BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintLinearDerivativeNonlocalTotalStrain1D)
#endif // ENABLE_SERIALIZATION
