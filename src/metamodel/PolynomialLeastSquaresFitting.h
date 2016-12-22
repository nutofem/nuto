#pragma once

// parent
#include "metamodel/Metamodel.h"


#include <vector>
#include "math/FullVector_Def.h"

namespace NuTo
{
class PolynomialLeastSquaresFitting : public virtual NuTo::Metamodel
{

public:
    //! @brief default constructor
    PolynomialLeastSquaresFitting();

    //! @brief ... Adds a pair of x and y coordinates that should be matched by the polynom
    //! @param rBoundaryCondition ... a pair of x and y coordinates that should be matched by the polynom
    void AddBoundaryCondition(std::pair<double,double> rBoundaryCondition);

    //! @brief ... Adds a pair of x and y coordinates that should be matched by the polynom
    //! @param rX ... the x value of the boundary condition
    //! @param rY ... the y value of the boundary condition
    void AddBoundaryCondition(double rX, double rY);

    //! @brief determine regression parameters
    void BuildDerived() override;

    //! @brief ... Gets the degree of the polynom
    //! @return degree of the polynom
    double GetDegree() const;

    //! @brief ... Gets the calculated polynomial coefficients
    //! @return ... polynomial coefficients
    NuTo::FullVector<double, Eigen::Dynamic> GetPolynomialCoefficients() const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId() const;

    //! @brief ... Sets the degree of the polynom
    //! @param rDegree ... degree of the polynom
    void SetDegree(int rDegree);

    //! @brief ... calculate approximation (in transformed space)
    //! @param rInputCoordinates ... matrix of input data points (transformed)
    //! @param rOutputCoordinates ... vector of output data (transformed)
    void SolveTransformed(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates)const;

protected:
    FullVector<double,Eigen::Dynamic>       mPolynomialCoeffs;
    std::vector<std::pair<double,double>>   mBoundaryConditions;
    int                                     mDegree             = -1;
};

}   // namespace NuTo

