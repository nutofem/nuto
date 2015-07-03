#ifndef STRUCTUREOUTPUTSPARSEMATRIXCSR
#define STRUCTUREOUTPUTSPARSEMATRIXCSR

#include "nuto/math/SparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBase.h"


namespace NuTo
{
//! @author Volker Hirthammer
//! @date June 25, 2015
//! @brief ...
class StructureOutputSparseMatrix: public StructureOutputBase
{
    friend class Structure;
public:
    StructureOutputSparseMatrix(SparseMatrix<double>& rMatrix)
        : StructureOutputBase(),
          mMatrix(rMatrix)
    {}

    virtual ~StructureOutputSparseMatrix(){}

    virtual int GetNumSubmatrices() const override
    {
        return 1;
    }

    virtual SparseMatrix<double>& GetSparseMatrixDouble() override
    {
        return mMatrix;
    }

private:
    bool mConstant  = false;
    bool mSymmetric = false;

    SparseMatrix<double>& mMatrix;
};


} // namespace NuTo


#endif // STRUCTUREOUTPUTSPARSEMATRIXCSR

