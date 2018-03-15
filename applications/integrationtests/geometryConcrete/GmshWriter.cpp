#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>
#include <cstdlib> // std::system
#include "geometryConcrete/GmshWriter.h"
#include "mechanics/mesh/MeshGmsh.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(Gmsh2D)
{
    boost::filesystem::path outputDir = boost::unit_test::framework::master_test_suite().argv[0] + std::string("Out");
    boost::filesystem::create_directory(outputDir);

    std::string gmshPath = boost::unit_test::framework::master_test_suite().argv[1];
    Eigen::MatrixXd m(2, 3);
    m.row(0) = Eigen::Vector3d(5, 5, 4);
    m.row(1) = Eigen::Vector3d(15, 15, 4);

    std::string filename = (outputDir / "2DInterface").string();

    auto opts = GmshWriter::Options::Default(2, 2, 0.5);
    opts.additionalOptions = "Mesh.SecondOrderIncomplete = 1;\n";

    GmshWriter::Write(filename + ".geo", GmshWriter::Box2D{Eigen::Vector2d(0, 0), Eigen::Vector2d(20, 20)}, m, opts);
    std::string cmd = gmshPath + " -2 -order 2 " + filename + ".geo -v 2";
    std::system(cmd.c_str());

    MeshGmsh meshGmsh(filename + ".msh");
    BOOST_CHECK(not meshGmsh.GetPhysicalGroup(opts.physicalGroupMatrix).Empty());
    BOOST_CHECK(not meshGmsh.GetPhysicalGroup(opts.physicalGroupAggregates).Empty());
    BOOST_CHECK(not meshGmsh.GetPhysicalGroup(opts.physicalGroupInterfaces).Empty());
}
