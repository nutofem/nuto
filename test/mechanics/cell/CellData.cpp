#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellData.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{
class Interpolation;
}


BOOST_AUTO_TEST_CASE(CellDataNodeValues)
{
    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
    fakeit::Mock<NuTo::Interpolation> interpolation0;
    NuTo::ElementSimple e0({&n0, &n1, &n2}, interpolation0.get());
    NuTo::DofType d0("dof0", 2, 0);

    NuTo::NodeSimple n3 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(1));
    NuTo::NodeSimple n4 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(3));
    NuTo::NodeSimple n5 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(7));
    fakeit::Mock<NuTo::Interpolation> interpolation1;
    NuTo::ElementSimple e1({&n3, &n4, &n5}, interpolation1.get());
    NuTo::DofType d1("dof0", 1, 1);

    NuTo::DofContainer<NuTo::ElementSimple*> elements;
    elements[d0] = &e0;
    elements[d1] = &e1;

    NuTo::CellData cell(elements);
    NuTo::NodeValues nodeValues0 = cell.GetNodeValues(d0);
    NuTo::NodeValues nodeValues1 = cell.GetNodeValues(d1);

    BoostUnitTest::CheckVector(nodeValues0, std::vector<double>({1,1,5,1,1,7}), 6);
    BoostUnitTest::CheckVector(nodeValues1, std::vector<double>({1,3,7}), 3);
}
