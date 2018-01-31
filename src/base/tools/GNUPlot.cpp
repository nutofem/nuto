#include "base/tools/GNUPlot.h"

#include "base/Exception.h"
#include "base/tools/PlotData.h"

#include <iostream>

#define GNUPLOTDATAFILENAME std::string("tmpGNUPlotData")

NuTo::Plot::GNUPlot::GNUPlot(bool keepOpen)
{
    if (keepOpen)
        mStream = popen("gnuplot -persistent", "w");
    else
        mStream = popen("gnuplot", "w");
}

NuTo::Plot::GNUPlot::~GNUPlot()
{
    if (mStream)
        fclose(mStream);
}


void NuTo::Plot::GNUPlot::AddPlot(Eigen::VectorXd x, Eigen::VectorXd y, std::array<unsigned char, 3> lineColor,
                                  eLineType lineType, std::string title)
{
    assert(x.rows() == y.rows());
    assert(x.cols() == y.cols());
    assert(x.cols() == 1);

    std::string tmpFileName = GNUPLOTDATAFILENAME + std::to_string(mSetupCurrentPlots.size()) + ".dat";

    std::FILE* dataFile = fopen(tmpFileName.c_str(), "w");
    if (!dataFile)
        throw Exception(__PRETTY_FUNCTION__, "Could not open temporary data file!");

    mSetupCurrentPlots.push_back({title, lineColor, lineType});
    std::string DataString;
    for (unsigned int i = 0; i < x.rows(); ++i)
    {
        DataString.append(std::to_string(x[i]));
        DataString.append(" ");
        DataString.append(std::to_string(y[i]));
        DataString.append("\n");
    }
    fwrite(DataString.c_str(), sizeof(char), DataString.length(), dataFile);

    fclose(dataFile);
}

void NuTo::Plot::GNUPlot::AddPlot(std::vector<double> x, std::vector<double> y, std::array<unsigned char, 3> lineColor,
                                  NuTo::Plot::eLineType lineType, std::string title)
{
    assert(x.size() == y.size());

    std::string tmpFileName = GNUPLOTDATAFILENAME + std::to_string(mSetupCurrentPlots.size()) + ".dat";

    std::FILE* dataFile = fopen(tmpFileName.c_str(), "w");
    if (!dataFile)
        throw Exception(__PRETTY_FUNCTION__, "Could not open temporary data file!");

    mSetupCurrentPlots.push_back({title, lineColor, lineType});
    std::string DataString;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        DataString.append(std::to_string(x[i]));
        DataString.append(" ");
        DataString.append(std::to_string(y[i]));
        DataString.append("\n");
    }
    fwrite(DataString.c_str(), sizeof(char), DataString.length(), dataFile);

    fclose(dataFile);
}


void NuTo::Plot::GNUPlot::Show() const
{
    if (mStream)
    {
        //        std::string test = "plot 'tmpGNUPlotData.dat' w linespoints \n";
        std::string command = "plot";
        for (unsigned int i = 0; i < mSetupCurrentPlots.size(); ++i)
        {
            command.append(" '")
                    .append(GNUPLOTDATAFILENAME)
                    .append(std::to_string(i))
                    .append(".dat' ")
                    .append(mSetupCurrentPlots[i].GetSetupString());
            if (i < mSetupCurrentPlots.size() - 1)
                command.append(", ");
        }
        command.append(" \n");

        fwrite(command.c_str(), sizeof(char), command.length(), mStream);
        fflush(mStream);
    }
}

void NuTo::Plot::GNUPlot::Clear()
{
    mSetupCurrentPlots.clear();
}
