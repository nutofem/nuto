#pragma once

#include <vector>
#include <string>
#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @brief Collection of helper functions for Eigen matrices
class EigenCompanion
{
public:
    //! @brief Append the rows of one matrix to the other
    //! @param top Matrix that is getting extended
    //! @param bottom Matrix that is added to the bottom of the other one
    static void AppendRows(Eigen::MatrixXd& top, const Eigen::MatrixXd& bottom);


    //! @brief writes a matrix to a file
    //! @param rMatrix matrix
    //! @param rFileName file name
    //! @param rDelimiter delimiters between the entries in each line, default = space
    static void WriteToFile(const Eigen::MatrixXd& rMatrix, const std::string& rFileName, std::string rDelimiter = " ");


    //! @brief reads a matrix from a file
    //! @param rFileName file name
    //! @param rDelimiter delimiters between the entries in each line, default = space
    //! @return matrix
    static Eigen::MatrixXd ReadFromFile(const std::string& rFileName);

    //! @brief converts data to a 3D vector, fills with zeros if needed
    //! @param data vector of arbitrary size
    //! @return 3D vector
    static Eigen::Vector3d To3D(const Eigen::VectorXd& data);

private:
    //! @brief converts a string, seperated by rDelimiter, to a vector of doubles
    //! @param rString string to convert
    //! @param rDelimiter delimiters between the entries in each line
    //! @return vector containing the numbers in the string
    static std::vector<double> StringToDoubles(const std::string& rString);
};
}
