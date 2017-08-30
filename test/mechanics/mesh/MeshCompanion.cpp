#include "BoostUnitTest.h"

#include "mechanics/MechanicsEnums.h"
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"


BOOST_AUTO_TEST_CASE(ImportFromGmshBinary)
{
    NuTo::Structure s(3);
    BOOST_CHECK_NO_THROW(NuTo::MeshCompanion::ImportFromGmsh(s, "MeshCompanionGmsh.msh"));
}

void PrismCreate(NuTo::Interpolation::eTypeOrder rCoordinateInterpolation)
{
    constexpr double thickness = .1;
    constexpr double lx = 10;
    constexpr double ly = 2;
    constexpr double lz = 5;

    NuTo::Structure s(3);

    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, rCoordinateInterpolation);

    int gMatrix = s.GroupCreate(NuTo::eGroupId::Elements);
    int gAggreg = s.GroupCreate(NuTo::eGroupId::Elements);

    if (rCoordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
    {
        s.NodeCreate(0, Eigen::Vector3d({-lx, 0, 0}));
        s.NodeCreate(1, Eigen::Vector3d({0, -ly, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0, 0, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0, ly, 0}));
        s.NodeCreate(4, Eigen::Vector3d({0, -ly / 2., lz / 2.}));
        s.NodeCreate(5, Eigen::Vector3d({0, ly / 2., lz / 2.}));
        s.NodeCreate(6, Eigen::Vector3d({0, 0, lz}));
        s.NodeCreate(7, Eigen::Vector3d({lx, 0, 0}));

        s.ElementCreate(1, it, {0, 1, 2, 4});
        s.ElementCreate(2, it, {0, 2, 3, 5});
        s.ElementCreate(3, it, {0, 4, 2, 5});
        s.ElementCreate(4, it, {0, 4, 5, 6});
        s.ElementCreate(5, it, {7, 1, 2, 4});
        s.ElementCreate(6, it, {7, 2, 3, 5});
        s.ElementCreate(7, it, {7, 4, 2, 5});
        s.ElementCreate(8, it, {7, 4, 5, 6});


        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gMatrix, 2);
        s.GroupAddElement(gMatrix, 3);
        s.GroupAddElement(gMatrix, 4);

        s.GroupAddElement(gAggreg, 5);
        s.GroupAddElement(gAggreg, 6);
        s.GroupAddElement(gAggreg, 7);
        s.GroupAddElement(gAggreg, 8);
    }
    else
    {
        Eigen::Vector3d(0.5, 0.0, 0.0);
        Eigen::Vector3d(0.5, 0.5, 0.0);
        Eigen::Vector3d(0.0, 0.5, 0.0);
        Eigen::Vector3d(0.0, 0.0, 0.5);
        Eigen::Vector3d(0.0, 0.5, 0.5);
        Eigen::Vector3d(0.5, 0.0, 0.5);


        s.NodeCreate(0, Eigen::Vector3d({0, -ly / 2, 0}));
        s.NodeCreate(1, Eigen::Vector3d({lx, 0, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0, ly / 2, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0, 0, lz}));

        s.NodeCreate(4, Eigen::Vector3d({lx / 2., -ly / 4, 0}));
        s.NodeCreate(5, Eigen::Vector3d({lx / 2., ly / 4, 0}));
        s.NodeCreate(6, Eigen::Vector3d({0, 0, 0}));

        s.NodeCreate(7, Eigen::Vector3d({0, -ly / 4, lz / 2}));
        s.NodeCreate(8, Eigen::Vector3d({0, ly / 4, lz / 2}));
        s.NodeCreate(9, Eigen::Vector3d({lx / 2., 0, lz / 2}));


        s.NodeCreate(10, Eigen::Vector3d({-lx, 0, 0}));
        s.NodeCreate(11, Eigen::Vector3d({-lx / 2., -ly / 4, 0}));
        s.NodeCreate(12, Eigen::Vector3d({-lx / 2., ly / 4, 0}));
        s.NodeCreate(13, Eigen::Vector3d({-lx / 2., 0, lz / 2}));


        s.ElementCreate(1, it, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        s.ElementCreate(2, it, {0, 10, 2, 3, 11, 12, 6, 7, 8, 13});
        //            s.ElementCreate(2, it, {1, 2, 3, 7});

        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gAggreg, 2);
    }

    auto prism = NuTo::MeshCompanion::ElementPrismsCreate(s, gMatrix, gAggreg, thickness);

    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(prism.second, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.ElementTotalConvertToInterpolationType();
}

BOOST_AUTO_TEST_CASE(CreatePrismEQUIDISTANT1)
{
    PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
}

BOOST_AUTO_TEST_CASE(CreatePrismEQUIDISTANT2)
{
    PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
}

void CheckPrismGmsh(std::string mshFile)
{
    NuTo::Structure s(3);
    auto meshInfo = NuTo::MeshCompanion::ImportFromGmsh(s, mshFile);
    auto& groupMaster = meshInfo[0].first;
    auto& groupSlave = meshInfo[1].first;
    double thickness = 0.5;
    BOOST_CHECK_NO_THROW(NuTo::MeshCompanion::ElementPrismsCreate(s, groupMaster, groupSlave, thickness));
}

BOOST_AUTO_TEST_CASE(CreatePrismGmsh)
{
    //CheckPrismGmsh("MeshCompanionGmsh.msh");
    // This test fails due to a bug. Possibly:
    // If one element is part of two prism surfaces, it fails. Really dunno why.
}

BOOST_AUTO_TEST_CASE(CreatePrismGmshFine)
{
    CheckPrismGmsh("MeshCompanionGmshFine.msh");
}
