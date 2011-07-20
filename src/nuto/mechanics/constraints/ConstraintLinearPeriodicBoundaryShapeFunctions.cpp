// $Id: ConstraintGlobalCrackOpening.cpp 314 2010-09-27 16:31:43Z unger3 $

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
#include "nuto/mechanics/constraints/ConstraintLinearPeriodicBoundaryShapeFunctions.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"

// constructor
NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::ConstraintLinearPeriodicBoundaryShapeFunctions(const StructureMultiscale* rStructure, int rPeriodicShapeFunction, double rValue):
        ConstraintLinear()
{
    mStructure = rStructure;
    if (rPeriodicShapeFunction<0 || rPeriodicShapeFunction>2)
    	throw MechanicsException("[NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::ConstraintLinearPeriodicBoundaryShapeFunctions] there are only 3 shape functions (0,1,2).");
    mPeriodicShapeFunction = rPeriodicShapeFunction;
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    rRHS(curConstraintEquation,0) = mRHS;

    rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDOFPeriodicBoundaryDisplacements()[mPeriodicShapeFunction],1);

    curConstraintEquation++;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearPeriodicBoundaryShapeFunctions" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mPeriodicShapeFunction)
       & BOOST_SERIALIZATION_NVP(mRHS);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearPeriodicBoundaryShapeFunctions" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions)
#endif // ENABLE_SERIALIZATION
