//
// Created by Thomas Titscher on 10/28/16.
//
#pragma once
#include <string>
#include <eigen3/Eigen/Core>
#include "mechanics/MechanicsException.h"
#include "../../math/Interpolation.h"
#include "math/EigenCompanion.h"

namespace NuTo
{
namespace Tools
{

//! @brief class that helps with the evaluation of global load displacement curves
class GlobalFractureEnergyIntegrator
{
public:

    //! @brief constructs from data vectors
    //! @param rForces vector containing forces
    //! @param rDisplacements vector container displacements
    //! @return obj
    GlobalFractureEnergyIntegrator(const Eigen::VectorXd& rForces, const Eigen::VectorXd& rDisplacements)
    {
        mForce = rForces.cwiseAbs();
        mDispl = rDisplacements.cwiseAbs();
        CheckData();
    }

    //! @brief constructs from data files
    //! @param rForces <path>/<file>.dat to container with force data
    //! @param rDisplacements <path>/<file>.dat to container with force data
    //! @param rColumn column to read from, hint for most cases: 1=x, 2=y, 3=z direction
    //! @return obj
    GlobalFractureEnergyIntegrator(std::string rFileForces, std::string rFileDisplacements, int rColumn)
    {
        mForce = NuTo::EigenCompanion::ReadFromFile(rFileForces).col(rColumn).cwiseAbs();
        mDispl = NuTo::EigenCompanion::ReadFromFile(rFileDisplacements).col(rColumn).cwiseAbs();

        CheckData();
    }


    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rArea "crack" area
    //! @return global GF
    double IntegrateSofteningCurve(double rArea)
    {
        return IntegrateSofteningCurveInternal(static_cast<int>(mForce.rows())) / rArea;
    }

    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rArea "crack" area
    //! @param rForceThreshold integration stops if the force drops below rForceThreshold
    //! @return global GF
    double IntegrateSofteningCurve(double rArea, double rForceThreshold)
    {
        int iEnd = FindIndexJustAboveForceThreshold(rForceThreshold);
        double integral = IntegrateSofteningCurveInternal(iEnd);
        integral += IntegrateIndexToForceThreshold(iEnd, rForceThreshold);
        return integral / rArea;
    }

private:

    //! @brief integrates F[rIndexEnd]
    //! @param rIndexEnd
    //! @param rForceThreshold
    //! @return
    double IntegrateIndexToForceThreshold(int rIndexEnd, double rForceThreshold)
    {
        double d0 = mDispl(rIndexEnd - 1);
        double d1 = mDispl(rIndexEnd);
        double f0 = mForce(rIndexEnd - 1);
        double f1 = mForce(rIndexEnd);

        double dX = d0 - d1;
        double dF = f0 - f1;
        double dThreshold  = d0 + dX / dF * (rForceThreshold - f0);
        //ToDo:: This interpolation should be done with NuTo::Math::Interpolation but I honestly have no idea how that class works. TT.

        double meanForceStop = (rForceThreshold + f0) / 2.;
        double deltaDisplStop = dThreshold - d0;

        return meanForceStop * deltaDisplStop;
    }

    //! @brief finds the index where the force is just before the rForceThreshold
    //! @param rForceThreshold force threshold
    //! @return index of F[index] > rForceThreshold && F[index+1] < rForceThreshold
    int FindIndexJustAboveForceThreshold(double rForceThreshold)
    {
        auto numRows = mForce.rows();
        if (rForceThreshold < mForce[numRows - 1])
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "The forces did not drop below the force threshold");

        for (auto i = 0; i < numRows-1; ++i)
        {
            if (mForce[i] > rForceThreshold && mForce[i+1] <= rForceThreshold)
                return i;
        }
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "The forces never reach the force threshold");
    }

    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rEnd integrate from index=0 to index=rEnd
    //! @return global GF
    double IntegrateSofteningCurveInternal(int rEnd)
    {
        double integral = 0.;
        for (int iRow = 0; iRow < rEnd - 1; ++iRow)
        {
            // integrate normally
            double meanForce = (mForce(iRow + 1) + mForce(iRow)) / 2.;
            double deltaDispl = mDispl(iRow + 1) - mDispl(iRow);
            integral += meanForce * deltaDispl;
        }
        return integral;
    }


    void CheckData() const
    {
        if (mForce.rows() != mDispl.rows())
            throw MechanicsException(__PRETTY_FUNCTION__, "Force.Rows != Displacement.Rows()");

        if (std::abs(mForce[0]) > 1.e-10)
            throw MechanicsException(__PRETTY_FUNCTION__, "Forces[0] != 0 ( but has to.)");

        if (std::abs(mDispl[0]) > 1.e-10)
            throw MechanicsException(__PRETTY_FUNCTION__, "Displacements[0] != 0 ( but has to.)");

    }

    Eigen::VectorXd mForce;
    Eigen::VectorXd mDispl;
};

}   // namespace Tools
}   // namespace NuTo