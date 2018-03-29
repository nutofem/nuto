#include "BoostUnitTest.h"
#include "nuto/mechanics/cell/Assembler.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(VectorAssembly)
{
    NuTo::DofType d("d", 1);
    DofVector<double> v;
    v[d] = Eigen::Vector3d(1, 2, 3);

    DofVector<int> numbering0, numbering1;
    numbering0[d] = Eigen::Vector3i(0, 1, 2);
    numbering1[d] = Eigen::Vector3i(2, 3, 4);

    DofContainer<int> size;
    size[d] = 5;
    VectorAssembler vectorAssembler(size);
    vectorAssembler.Add(v, numbering0);
    vectorAssembler.Add(v, numbering1);

    BoostUnitTest::CheckVector(vectorAssembler.Get()[d], std::vector<double>{1, 2, 3 + 1, 2, 3}, 5);

    vectorAssembler.Reset();
    BoostUnitTest::CheckVector(vectorAssembler.Get()[d], std::vector<double>{0, 0, 0, 0, 0}, 5);
}
