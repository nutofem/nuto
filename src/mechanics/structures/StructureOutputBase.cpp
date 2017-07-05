
#include "base/Exception.h"
#include "StructureOutputBase.h"

NuTo::StructureOutputBase::StructureOutputBase()
{
}

NuTo::StructureOutputBase::~StructureOutputBase()
{
}

NuTo::StructureOutputBlockMatrix& NuTo::StructureOutputBase::AsStructureOutputBlockMatrix()
{
    throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "[ StructureOutput is not of type BlockMatrix" );
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBase::AsStructureOutputBlockVector()
{
    throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "[ StructureOutput is not of type BlockVector" );
}

void NuTo::StructureOutputBase::SetSymmetry(bool rSymmetric)
{
    throw Exception("[StructureOutputBase::SetSymmetry] symmetry is not stored.");
}

bool NuTo::StructureOutputBase::IsSymmetric() const
{
    throw Exception("[StructureOutputBase::SetSymmetry] symmetry is not stored.");
}

void NuTo::StructureOutputBase::SetZero()
{
    throw Exception("[StructureOutputBase::SetZero] not implemented.");
}
