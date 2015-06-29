#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "StructureOutputBase.h"

NuTo::StructureOutputBase::StructureOutputBase()
{}

NuTo::StructureOutputBase::~StructureOutputBase()
{}

std::shared_ptr<NuTo::SparseMatrixCSRVector2<double> > &NuTo::StructureOutputBase::GetPtrSparseMatrixCSRVector2Double()
{
    throw MechanicsException("[StructureOutputBase::GetPtrSparseMatrixCSRVector2Double] output matrix is not of type std::shared_ptr<SparseMatrixCSRVector2<double> >");
}

NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& NuTo::StructureOutputBase::GetFullMatrixDouble()
{
    throw MechanicsException("[StructureOutputBase::GetFullMatrixDouble] output matrix is not of type FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>");
}

NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& NuTo::StructureOutputBase::GetFullMatrixInt()
{
    throw MechanicsException("[StructureOutputBase::GetFullMatrixDouble] output matrix is not of type FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>");
}

NuTo::FullVector<double,Eigen::Dynamic>& NuTo::StructureOutputBase::GetFullVectorDouble()
{
    throw MechanicsException("[StructureOutputBase::GetFullVectorDouble] output matrix is not of type FullVector<double,Eigen::Dynamic>");
}

NuTo::FullVector<int,Eigen::Dynamic>& NuTo::StructureOutputBase::GetFullVectorInt()
{
    throw MechanicsException("[StructureOutputBase::GetFullVectorInt] output matrix is not of type FullVector<double,Eigen::Dynamic>");
}

std::vector<int>& NuTo::StructureOutputBase::GetVectorInt()
{
    throw MechanicsException("[StructureOutputBase::GetFullMatrixDouble] output matrix is not of type std::vector<int>");
}

void NuTo::StructureOutputBase::SetSymmetry(bool rSymmetric)
{
    throw MechanicsException("[StructureOutputBase::SetSymmetry] symmetry is not stored.");
}

bool NuTo::StructureOutputBase::GetSymmetry()const
{
    throw MechanicsException("[StructureOutputBase::SetSymmetry] symmetry is not stored.");
}

void NuTo::StructureOutputBase::SetConstant(bool rConstant)
{
    throw MechanicsException("[StructureOutputBase::SetConstant] constness is not stored.");
}

bool NuTo::StructureOutputBase::GetConstant()const
{
    throw MechanicsException("[StructureOutputBase::GetConstant] constness is not stored.");
}
