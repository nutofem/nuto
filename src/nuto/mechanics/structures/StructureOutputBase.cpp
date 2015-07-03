#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "StructureOutputBase.h"

NuTo::StructureOutputBase::StructureOutputBase()
{}

NuTo::StructureOutputBase::~StructureOutputBase()
{}

NuTo::SparseMatrix<double> &NuTo::StructureOutputBase::GetSparseMatrixDouble()
{
    throw MechanicsException("[StructureOutputBase::GetSparseMatrixDouble] output matrix is not of type NuTo::SparseMatrix<double>& >");
}

NuTo::SparseMatrix<double> &NuTo::StructureOutputBase::GetSparseMatrixDouble(NuTo::StructureEnum::eSubMatrix rSubmatrixEnum)
{
    throw MechanicsException("[StructureOutputBase::GetSparseMatrixDouble] output matrix is not of type NuTo::SparseMatrix<double>& >");
}


NuTo::FullVector<double,Eigen::Dynamic>& NuTo::StructureOutputBase::GetFullVectorDouble()
{
    throw MechanicsException("[StructureOutputBase::GetFullVectorDouble] output matrix is not of type FullVector<double,Eigen::Dynamic>");
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
