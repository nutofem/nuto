#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>
#include <cstdlib> // std::system
#include "nuto/geometryConcrete/GmshWriter.h"
#include "nuto/mechanics/mesh/MeshGmsh.h"

using namespace NuTo;

void CheckMesh(std::string meshFile, GmshWriter::Options opts)
{
    MeshGmsh meshGmsh(meshFile);
    BOOST_CHECK(not meshGmsh.GetPhysicalGroup(opts.physicalGroupMatrix).Empty());
    BOOST_CHECK(not meshGmsh.GetPhysicalGroup(opts.physicalGroupAggregates).Empty());

    if (opts.interfaceThickness != 0)
        BOOST_CHECK(not meshGmsh.GetPhysicalGroup(opts.physicalGroupInterfaces).Empty());
}

std::string GmshPath()
{
    return boost::unit_test::framework::master_test_suite().argv[1];
}

std::string FileName(std::string testName)
{
    boost::filesystem::path outputDir = boost::unit_test::framework::master_test_suite().argv[0] + std::string("Out");
    boost::filesystem::create_directory(outputDir);
    return (outputDir / testName).string();
}

template <typename TBox, typename TAggregates>
void Check(std::string filename, TBox box, TAggregates aggregates, GmshWriter::Options opts)
{
    GmshWriter::Write(filename + ".geo", box, aggregates, opts);
    int dim = aggregates.cols() - 1;
    std::system((GmshPath() + " -" + std::to_string(dim) + " -order 2 " + filename + ".geo -v 2").c_str());
    // -v 2 supresses normal output
    CheckMesh(filename + ".msh", opts);
}

BOOST_AUTO_TEST_CASE(Gmsh2D)
{
    Eigen::MatrixXd m(2, 3);
    m.row(0) = Eigen::Vector3d(5, 5, 3);
    m.row(1) = Eigen::Vector3d(15, 15, 3);

    GmshWriter::Options opts(2, 2);
    Check(FileName("2DNormal"), GmshWriter::Box2D{Eigen::Vector2d(0, 0), Eigen::Vector2d(20, 20)}, m, opts);


    opts.interfaceThickness = 0.5;
    opts.additionalOptions = "Mesh.SecondOrderIncomplete = 1;\n";
    Check(FileName("2DInterface"), GmshWriter::Box2D(20, 20), m, opts);
}

BOOST_AUTO_TEST_CASE(Gmsh3D)
{
    Eigen::MatrixXd m(2, 4);
    m.row(0) = Eigen::Vector4d(5, 5, 5, 3);
    m.row(1) = Eigen::Vector4d(15, 15, 15, 3);

    GmshWriter::Options opts(1.25, 1.25);
    Check(FileName("3DBoxNormal"), GmshWriter::Box3D{Eigen::Vector3d::Zero(), Eigen::Vector3d::Constant(20)}, m, opts);


    opts.interfaceThickness = 0.5;
    Check(FileName("3DBoxInterface"), GmshWriter::Box3D{20, 20, 20}, m, opts);
}

BOOST_AUTO_TEST_CASE(Gmsh3DCylinder)
{
    Eigen::MatrixXd m(2, 4);
    m.row(0) = Eigen::Vector4d(0, 0, 5, 3);
    m.row(1) = Eigen::Vector4d(0, 0, -5, 3);

    GmshWriter::Options opts(1.3, 1.3);
    Check(FileName("3DCylinderNormal"), GmshWriter::Cylinder{5, 10}, m, opts);


    opts.interfaceThickness = 0.5;
    Check(FileName("3DCylinderInterface"), GmshWriter::Cylinder{5, 10}, m, opts);
}
