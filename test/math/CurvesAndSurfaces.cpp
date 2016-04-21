#include "nuto/math/BSplineCurve.h"
#include <boost/filesystem.hpp>
#include <math.h>

int main(int argc,char *argv[])
{
    // point sequence to interpolate
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> points;

    int numPoints = 10;
    int dim = 2;
    int degree = 3;

    points.Resize(numPoints, dim);

    // half a circle
    double phi = 0.;
    for (int i = 0; i < numPoints-1; i++)
    {
        points.SetValue(i,0,std::cos(phi));
        points.SetValue(i,1,std::sin(phi));
        phi+=M_PI/numPoints;
    }
    points.SetValue(numPoints-1,0,-1.);
    points.SetValue(numPoints-1,1,0.);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> AInv;

    NuTo::BSplineCurve curve(degree, points, AInv);

    int num = 100;
    NuTo::FullVector<double, Eigen::Dynamic> parameters(num);
    double inc = (1./num);
    for(int i = 1; i < num-1; i++) parameters[i] = i*inc;
    parameters[num-1]=1.;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> curvepoints = curve.CurvePoints(parameters);

    // write the data
    std::string resultDir = std::string("./CurvesAndSurfacesResult");
    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDir);

    boost::filesystem::path resultFileName(resultDir);
    std::string ident("curve");
    resultFileName /= ident+".dat";

    curvepoints.WriteToFile(resultFileName.string(), "  ");

    return 0;
}
