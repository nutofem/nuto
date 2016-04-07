#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "StructureOutputBase.h"

NuTo::StructureOutputBase::StructureOutputBase()
{}

NuTo::StructureOutputBase::~StructureOutputBase()
{}

NuTo::StructureOutputBlockMatrix& NuTo::StructureOutputBase::AsStructureOutputBlockMatrix()
{
    throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "[ StructureOutput is not of type BlockMatrix" );
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBase::AsStructureOutputBlockVector()
{
    throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "[ StructureOutput is not of type BlockVector" );
}

void NuTo::StructureOutputBase::SetSymmetry(bool rSymmetric)
{
    throw MechanicsException("[StructureOutputBase::SetSymmetry] symmetry is not stored.");
}

bool NuTo::StructureOutputBase::IsSymmetric()const
{
    throw MechanicsException("[StructureOutputBase::SetSymmetry] symmetry is not stored.");
}

void NuTo::StructureOutputBase::SetZero()
{
    throw MechanicsException("[StructureOutputBase::SetZero] not implemented.");
}
