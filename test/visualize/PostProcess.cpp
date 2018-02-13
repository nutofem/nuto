#include "visualize/AverageGeometries.h"
#include "visualize/PostProcess.h"
#include "TestStructure.h"

#include <boost/filesystem.hpp>

using namespace NuTo;
using namespace NuTo::Visualize;
namespace fs = boost::filesystem;

BOOST_AUTO_TEST_CASE(PostProcessResultDir)
{
    PostProcess postprocess("PostProcessOutput");
    BOOST_CHECK(fs::exists(postprocess.ResultDirectory()));
    BOOST_CHECK(fs::is_directory(postprocess.ResultDirectory()));

    postprocess.ResultDirectory("PostProcessOutput/");
    BOOST_CHECK(fs::exists(postprocess.ResultDirectory()));
    BOOST_CHECK(fs::is_directory(postprocess.ResultDirectory()));
}

BOOST_AUTO_TEST_CASE(PostProcessTest)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    PostProcess postprocess("PostProcessOutput/");
    postprocess.DefineVisualizer("blub", cells, AverageHandler(AverageGeometryQuad()));
    postprocess.Plot(42., false);
    postprocess.Plot(6174., false);

    fs::path expectedPvdFile = fs::path(postprocess.ResultDirectory()) / "blub.pvd";
    fs::path expectedvtuFile0 = fs::path(postprocess.ResultDirectory()) / "blub0.vtu";
    fs::path expectedvtuFile1 = fs::path(postprocess.ResultDirectory()) / "blub1.vtu";

    BOOST_CHECK(fs::is_regular_file(expectedPvdFile));
    BOOST_CHECK(fs::is_regular_file(expectedvtuFile0));
    BOOST_CHECK(fs::is_regular_file(expectedvtuFile1));

    pt::ptree tree;
    pt::read_xml(expectedPvdFile.string(), tree);
    double timeStep0 = tree.get<double>("VTKFile.Collection.DataSet.<xmlattr>.timestep");
    std::string file0 = tree.get<std::string>("VTKFile.Collection.DataSet.<xmlattr>.file");

    BOOST_CHECK_CLOSE(timeStep0, 42., 1.e-10);
    BOOST_CHECK(file0 == "blub0.vtu");
}

BOOST_AUTO_TEST_CASE(PostProcessEdgyMcEdge)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    // Add visualize object before defining the visualizer
    DofType dof("stuff", 1);
    PostProcess postprocess;
    BOOST_CHECK_THROW(postprocess.Add("edge", dof), Exception);

    postprocess.DefineVisualizer("edge", cells, AverageHandler(AverageGeometryQuad()));

    // defining the same visualizer twice
    BOOST_CHECK_THROW(postprocess.DefineVisualizer("edge", cells, AverageHandler(AverageGeometryQuad())), Exception);

    // call `Plot` without setting a result directory
    BOOST_CHECK_THROW(postprocess.Plot(0), Exception);

    // call `ResultDirectory` with an emptry string
    BOOST_CHECK_THROW(postprocess.ResultDirectory(""), Exception);
}
