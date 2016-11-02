//
// Created by Thomas Titscher on 10/28/16.
//
#pragma once
#include "nuto/math/FullVector.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "../../math/Interpolation.h"
#include <string>

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
    GlobalFractureEnergyIntegrator(const FullVector<double, Eigen::Dynamic>& rForces, const FullVector<double, Eigen::Dynamic>& rDisplacements)
        : mForce(rForces), mDispl(rDisplacements)
    {
        CheckData();
    }

    //! @brief constructs from data files
    //! @param rForces <path>/<file>.dat to container with force data
    //! @param rDisplacements <path>/<file>.dat to container with force data
    //! @param rColumn column to read from, hint for most cases: 1=x, 2=y, 3=z direction
    //! @return obj
    GlobalFractureEnergyIntegrator(std::string rFileForces, std::string rFileDisplacements, int rColumn)
    {
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> reader;

        reader.ReadFromFile(rFileForces);
        mForce = reader.GetColumn(rColumn);

        reader.ReadFromFile(rFileDisplacements);
        mDispl = reader.GetColumn(rColumn);

        CheckData();
    }


    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rArea "crack" area
    //! @param rForceThreshold integration stops if the force drops below rForceThreshold
    //! @return global GF
    double IntegrateSofteningCurve(double rArea, double rForceThreshold)
    {
        unsigned numRows = mForce.rows();
        double integral = 0.;
        double fLast = mForce[numRows - 1];

        if (rForceThreshold < fLast)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "The forces did not drop below the force threshold");

        bool isPeakReached = false;
        unsigned int iRow = 1;
        while (not isPeakReached or mForce(iRow) > rForceThreshold)  // iRow < 10 (or so...) to fully include the elastic part
        {
            // integrate normally
            double meanForce = (mForce(iRow) + mForce(iRow - 1)) / 2.;
            double deltaDispl = mDispl(iRow) - mDispl(iRow - 1);

            integral += meanForce * deltaDispl;
            ++iRow;
            if (mForce(iRow) - mForce(iRow - 1) < 0)
                isPeakReached = true;

            if (iRow - 1 == numRows)
                throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "The forces never reach the force threshold");
        }

        // force(iRow) just dropped below rForceThreshold
        // now integrate the missing part from the last mForce until rForceThreshold
        double d0 = mDispl(iRow - 1);
        double d1 = mDispl(iRow);
        double f0 = mForce(iRow - 1);
        double f1 = mForce(iRow);

        if (f0 < rForceThreshold) // just to be 100% sure
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "f0 < force threshold");

        if (f1 > rForceThreshold) // just to be 100% sure
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "f1 > force threshold");


        double dX = d0 - d1;
        double dF = f0 - f1;
        double dThreshold  = d0 + dX / dF * (rForceThreshold - f0);
        //ToDo:: This interpolation should be done with NuTo::Math::Interpolation but I honestly have no idea how that class works. TT.

        double meanForceStop = (rForceThreshold + f0) / 2.;
        double deltaDisplStop = dThreshold - d0;

        integral += meanForceStop * deltaDisplStop;

        return integral / rArea;
    }

private:

    void CheckData() const
    {
        if (mForce.rows() != mDispl.rows())
            throw MechanicsException(__PRETTY_FUNCTION__, "Force.Rows != Displacement.Rows()");

        if (std::abs(mForce[0]) > 1.e-10)
            throw MechanicsException(__PRETTY_FUNCTION__, "Forces[0] != 0 ( but has to.)");

        if (std::abs(mDispl[0]) > 1.e-10)
            throw MechanicsException(__PRETTY_FUNCTION__, "Displacements[0] != 0 ( but has to.)");

    }

    FullVector<double, Eigen::Dynamic> mForce;
    FullVector<double, Eigen::Dynamic> mDispl;
};

}   // namespace Tools
}   // namespace NuTo