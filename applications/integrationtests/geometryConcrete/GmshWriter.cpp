#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>
#include <cstdlib> // std::system
#include "geometryConcrete/GmshWriter.h"
#include "mechanics/mesh/MeshGmsh.h"

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
    std::system((GmshPath() + " -2 -order 2 " + filename + ".geo -v 2").c_str());
    CheckMesh(filename + ".msh", opts);
}

BOOST_AUTO_TEST_CASE(Gmsh2D)
{
    Eigen::MatrixXd m(2, 3);
    m.row(0) = Eigen::Vector3d(5, 5, 4);
    m.row(1) = Eigen::Vector3d(15, 15, 4);

    auto opts = GmshWriter::Options::Default(2, 2, 0.5);
    opts.additionalOptions = "Mesh.SecondOrderIncomplete = 1;\n";

    Check(FileName("2DInterface"), GmshWriter::Box2D{Eigen::Vector2d(0, 0), Eigen::Vector2d(20, 20)}, m, opts);
}
