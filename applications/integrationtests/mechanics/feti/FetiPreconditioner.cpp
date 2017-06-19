
#include "mechanics/feti/FetiIdentityPreconditioner.h"
#include "mechanics/MechanicsException.h"
#include <boost/test/unit_test.hpp>

using VectorType = NuTo::FetiPreconditioner::VectorType;

constexpr int vectorSize = 12345;
const VectorType vec = VectorType::Random(vectorSize);

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    NuTo::FetiIdentityPreconditioner preconditioner;

    if (not vec.isApprox(preconditioner.ApplyOnTheLeft(vec)))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Identity preconditioner changes input when applied.");

    MPI_Finalize();
}
