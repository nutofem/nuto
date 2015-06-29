#ifndef STRUCTUREOUTPUTSPARSEMATRIXCSR
#define STRUCTUREOUTPUTSPARSEMATRIXCSR

#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/mechanics/structures/StructureOutputBase.h"


namespace NuTo
{
//! @author Volker Hirthammer
//! @date June 25, 2015
//! @brief ...
class StructureOutputSparseMatrixCSRVector2 : public StructureOutputBase
{
    friend class Structure;
public:
    StructureOutputSparseMatrixCSRVector2() : StructureOutputBase(){}

    virtual ~StructureOutputSparseMatrixCSRVector2(){}

    virtual std::shared_ptr<SparseMatrixCSRVector2<double>>& GetPtrSparseMatrixCSRVector2Double() override
    {
        return mMatrix;
    }

private:
    bool mConstant  = false;
    bool mSymmetric = false;

    std::shared_ptr<SparseMatrixCSRVector2<double>> mMatrix;
};


} // namespace NuTo


#endif // STRUCTUREOUTPUTSPARSEMATRIXCSR

