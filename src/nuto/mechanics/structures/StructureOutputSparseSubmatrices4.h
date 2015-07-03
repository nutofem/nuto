#ifndef STRUCTUREOUTPUTSPARSESUBMATRICES4
#define STRUCTUREOUTPUTSPARSESUBMATRICES4


#include "nuto/math/SparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBase.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"


namespace NuTo
{
//! @author Volker Hirthammer
//! @date July 6, 2015
//! @brief ...
class StructureOutputSparseSubmatrices4: public StructureOutputBase
{
    friend class Structure;
public:
    StructureOutputSparseSubmatrices4(SparseMatrix<double>& rMatrixJJ,
                                      SparseMatrix<double>& rMatrixJK,
                                      SparseMatrix<double>& rMatrixKJ,
                                      SparseMatrix<double>& rMatrixKK)
        : StructureOutputBase(),
          mMatrixJJ(rMatrixJJ),
          mMatrixJK(rMatrixJK),
          mMatrixKJ(rMatrixKJ),
          mMatrixKK(rMatrixKK)
    {}

    virtual ~StructureOutputSparseMatrix(){}

    virtual int GetNumSubmatrices() const override
    {
        return 4;
    }

    virtual SparseMatrix<double>& GetSparseMatrixDouble(StructureEnum::eSubMatrix rSubmatrixEnum) override
    {
        switch(rSubmatrixEnum)
        {
            case StructureEnum::eSubMatrix::JJ:
            {
                return mMatrixJJ;
            }
            case StructureEnum::eSubMatrix::JK:
            {
                return mMatrixJK;
            }
            case StructureEnum::eSubMatrix::JJ:
            {
                return mMatrixKJ;
            }
            case StructureEnum::eSubMatrix::KJ:
            {
                return mMatrixKK;
            }
            default:
            {
                throw MechanicsException("[StructureOutputSparseSubmatrices4::GetSparseMatrixDouble] requested submatrix not implemented");
            }
        }
    }

private:
    bool mConstant  = false;
    bool mSymmetric = false;

    SparseMatrix<double>& mMatrixJJ, mMatrixJK, mMatrixKJ, mMatrixKK;
};


} // namespace NuTo

#endif // STRUCTUREOUTPUTSPARSESUBMATRICES4

