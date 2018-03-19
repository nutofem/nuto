#pragma once

#include <string>
#include <Eigen/Core>

namespace NuTo
{

//! @brief Collection of helper functions for Eigen matrices
class EigenIO
{
public:
    //! @brief writes a matrix to a file
    //! @param rMatrix matrix
    //! @param rFileName file name
    //! @param rDelimiter delimiters between the entries in each line, default = space
    static void WriteToFile(const Eigen::MatrixXd& rMatrix, const std::string& rFileName, std::string rDelimiter = " ");

    //! @brief reads a matrix from a file
    //! @param rFileName file name
    //! @return matrix
    static Eigen::MatrixXd ReadFromFile(const std::string& rFileName);
};
} /* NuTo */
