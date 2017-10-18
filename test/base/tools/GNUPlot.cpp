#include "base/tools/GNUPlot.h"


using namespace NuTo::Plot;

int main(int argc, char* argv[])
{
    GNUPlot gnuplot(true);

    Eigen::VectorXd x(5);
    Eigen::VectorXd y(5);

    x << 0, 1, 2, 3, 4;
    y << 0, 1, 2, 1, 0;

    gnuplot.AddPlot(x, y);

    y << 0, 3.5, 2.1, 1.5, 3.;

    gnuplot.AddPlot(x, y, {255, 0, 0}, eLineType::LINESPOINTS, "Data 2");
    gnuplot.Show();

    gnuplot.Clear();
    gnuplot.AddPlot(x, y, {255, 0, 255}, eLineType::LINESPOINTS, "Data 2");
    gnuplot.Show();
    return 0;
}
