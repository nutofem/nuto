#ifndef STRUCTUREOUTPUTBASE_H
#define STRUCTUREOUTPUTBASE_H

#include <nuto/math/SparseMatrixCSRVector2.h>
#include <nuto/mechanics/structures/StructureBaseEnum.h>

namespace NuTo
{
//! @author Volker Hirthammer
//! @date June 25, 2015
//! @brief ...
class StructureOutputBase
{
public:
    StructureOutputBase();

    virtual ~StructureOutputBase();

    virtual SparseMatrix<double>& GetSparseMatrixDouble();

    virtual SparseMatrix<double>& GetSparseMatrixDouble(StructureEnum::eSubMatrix rSubmatrixEnum);

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble();

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble(StructureEnum::eSubVector rSubvectorEnum);

    virtual int GetNumSubmatrices() const;

    virtual int GetNumSubvectors() const;

    virtual void SetSymmetry(bool rSymmetric);

    virtual bool GetSymmetry()const;

    virtual void SetConstant(bool rConstant);

    virtual bool GetConstant()const;
};

} // namespace NuTo

#endif // STRUCTUREOUTPUTBASE_H
