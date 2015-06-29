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

    virtual std::shared_ptr<SparseMatrixCSRVector2<double>>& GetPtrSparseMatrixCSRVector2Double();

    virtual FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixDouble();

    virtual FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixInt();

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble();

    virtual FullVector<int,Eigen::Dynamic>& GetFullVectorInt();

    virtual std::vector<int>& GetVectorInt();

    virtual void SetSymmetry(bool rSymmetric);

    virtual bool GetSymmetry()const;

    virtual void SetConstant(bool rConstant);

    virtual bool GetConstant()const;
};

} // namespace NuTo

#endif // STRUCTUREOUTPUTBASE_H
