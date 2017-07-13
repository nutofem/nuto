#include "math/EigenCompanion.h"

#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include "base/Exception.h"

using namespace NuTo;

void EigenCompanion::AppendRows(Eigen::MatrixXd& top, const Eigen::MatrixXd& bottom)
{
    if (top.cols() != bottom.cols())
        throw Exception(__PRETTY_FUNCTION__, "Number of columns for both matrices must be identical.");
    top.conservativeResize(top.rows() + bottom.rows(), Eigen::NoChange);
    top.bottomRows(bottom.rows()) = bottom;
}

void EigenCompanion::WriteToFile(const Eigen::MatrixXd& rMatrix, const std::string& rFileName, std::string rDelimiter)
{
    std::ofstream fileStream(rFileName.c_str());
    if (!fileStream.is_open())
        throw Exception(__PRETTY_FUNCTION__, "File " + rFileName + " could not be opened.");

    // go through the matrix and store the values
    fileStream.setf(std::ios::scientific, std::ios::floatfield);
    for (int iRow = 0; iRow < rMatrix.rows(); ++iRow)
    {
        for (int iCol = 0; iCol < rMatrix.cols(); ++iCol)
        {
            fileStream << std::setprecision(12) << std::setw(15) << std::showpoint << rMatrix(iRow, iCol);
            if (iCol != rMatrix.cols() - 1)
                fileStream << rDelimiter;
        }
        if (iRow != rMatrix.rows() - 1)
            fileStream << '\n';
    }
}

Eigen::MatrixXd EigenCompanion::ReadFromFile(const std::string& rFileName, std::string rDelimiter)
{
    std::ifstream fileStream(rFileName.c_str());
    if (!fileStream.is_open())
        throw Exception(__PRETTY_FUNCTION__, "File " + rFileName + " could not be opened.");

    // get total number of lines
    fileStream.seekg(0, std::ios::beg);
    unsigned int numRows = 0;
    std::string line;
    while (getline(fileStream, line))
        numRows++;

    // reset stream
    fileStream.clear();
    fileStream.seekg(0, std::ios::beg);

    // read first line in order to obtain number of columns
    getline(fileStream, line);
    unsigned int numColumns = StringToDoubles(line, rDelimiter).size();

    // resize matrix
    Eigen::MatrixXd m(numRows, numColumns);

    // reset stream
    fileStream.clear();
    fileStream.seekg(0, std::ios::beg);

    // read matrix
    for (unsigned int iRow = 0; iRow < numRows; ++iRow)
    {
        getline(fileStream, line);
        auto columnValues = StringToDoubles(line, rDelimiter);

        for (unsigned int iColumn = 0; iColumn < columnValues.size(); ++iColumn)
            m(iRow, iColumn) = columnValues[iColumn];
    }
    fileStream.close();
    return m;
}

Eigen::Vector3d EigenCompanion::To3D(const Eigen::VectorXd& data)
{
    const int dimension = data.rows();
    Eigen::Vector3d vector3d = Eigen::Vector3d::Zero();
    vector3d.block(0, 0, dimension, 1) = data;
    return vector3d;
}


std::vector<double> EigenCompanion::StringToDoubles(const std::string& rString, std::string rDelimiter)
{
    std::istringstream lineStream(rString);
    std::vector<std::string> splitLine;
    std::copy(std::istream_iterator<std::string>(lineStream), std::istream_iterator<std::string>(),
              std::back_inserter(splitLine));

    std::vector<double> doubles;
    doubles.reserve(splitLine.size());
    for (const auto& split : splitLine)
    {
        double d;
        std::istringstream(split) >> d;
        doubles.push_back(d); // std::stod throws [probably due to some result files that clash with the numeric_limits]
    }

    return doubles;
}
