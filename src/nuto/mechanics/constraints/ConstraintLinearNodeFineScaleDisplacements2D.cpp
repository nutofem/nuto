// $Id: ConstraintLinearNodeFineScaleDisplacements2D.cpp 309 2010-09-22 22:21:24Z unger3 $

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
#include "nuto/mechanics/nodes/NodeDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeFineScaleDisplacements2D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

NuTo::ConstraintLinearNodeFineScaleDisplacements2D::ConstraintLinearNodeFineScaleDisplacements2D(const NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue) :
        ConstraintNode(rNode), ConstraintLinear()
{
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeFineScaleDisplacements2D::ConstraintLinearNodeFineScaleDisplacements2D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.mEigenMatrix.data(),2*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeFineScaleDisplacements2D::ConstraintLinearNodeFineScaleDisplacements2D] direction vector has zero length.");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeFineScaleDisplacements2D::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix)const
{
    if (mNode->GetNumFineScaleDisplacements()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeFineScaleDisplacements2D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    if (fabs(mDirection[0])>1e-18)
        rConstraintMatrix.AddEntry(curConstraintEquation,mNode->GetDofFineScaleDisplacement(0),mDirection[0]);
    if (fabs(mDirection[1])>1e-18)
        rConstraintMatrix.AddEntry(curConstraintEquation,mNode->GetDofFineScaleDisplacement(1),mDirection[1]);

    curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::GetRHS(int& curConstraintEquation,NuTo::FullMatrix<double>& rRHS)const
{
    if (mNode->GetNumFineScaleDisplacements()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeFineScaleDisplacements2D::ConstraintBase] Node does not have displacements or has more than one displacement component.");

    rRHS(curConstraintEquation,0) = mRHS;

    curConstraintEquation++;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearNodeFineScaleDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeFineScaleDisplacements2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNode)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRHS)
       & BOOST_SERIALIZATION_NVP(mDirection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeFineScaleDisplacements2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeFineScaleDisplacements2D)
#endif // ENABLE_SERIALIZATION
