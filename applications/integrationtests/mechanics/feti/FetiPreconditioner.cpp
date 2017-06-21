
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/feti/FetiIdentityPreconditioner.h"
#include "mechanics/feti/FetiLumpedPreconditioner.h"
#include "mechanics/feti/FetiDirichletPreconditioner.h"
#include "mechanics/MechanicsException.h"

using VectorType = NuTo::FetiPreconditioner::VectorType;
using SparseMatrixType = NuTo::FetiPreconditioner::SparseMatrixType;

constexpr int vectorSize = 12345;
const VectorType vec = VectorType::Random(vectorSize);


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    SparseMatrixType B(1,1);

    NuTo::DofStatus dofStatus;
    NuTo::StructureOutputBlockMatrix hessian(dofStatus);

    NuTo::FetiIdentityPreconditioner identityPreconditioner;

    std::vector<int> lagrangeMultiplierDofIds;
    identityPreconditioner.Compute(hessian, B, lagrangeMultiplierDofIds);



    if (not vec.isApprox(identityPreconditioner.ApplyOnTheLeft(vec)))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Identity preconditioner changes input when applied.");

    NuTo::FetiLumpedPreconditioner lumpedPreconditioner;
    identityPreconditioner.Compute(hessian, B, lagrangeMultiplierDofIds);

    NuTo::FetiDirichletPreconditioner dirichletPreconditioner;
    identityPreconditioner.Compute(hessian, B, lagrangeMultiplierDofIds);


    MPI_Finalize();
}
