//
// Created by Thomas Titscher on 10/28/16.
//
#include "mechanics/tools/GlobalFractureEnergyIntegrator.h"
#include "base/Exception.h"
#include "math/EigenIO.h"

NuTo::Tools::GlobalFractureEnergyIntegrator::GlobalFractureEnergyIntegrator(const Eigen::VectorXd& rForces,
                                                                            const Eigen::VectorXd& rDisplacements)
{
    mForce = rForces.cwiseAbs();
    mDispl = rDisplacements.cwiseAbs();
    CheckData();
}

NuTo::Tools::GlobalFractureEnergyIntegrator::GlobalFractureEnergyIntegrator(std::string rFileForces,
                                                                            std::string rFileDisplacements, int rColumn)
{
    mForce = NuTo::EigenIO::ReadFromFile(rFileForces).col(rColumn).cwiseAbs();
    mDispl = NuTo::EigenIO::ReadFromFile(rFileDisplacements).col(rColumn).cwiseAbs();

    CheckData();
}


double NuTo::Tools::GlobalFractureEnergyIntegrator::IntegrateSofteningCurve(double rArea) const
{
    return IntegrateSofteningCurveInternal(static_cast<int>(mForce.rows())) / rArea;
}

double NuTo::Tools::GlobalFractureEnergyIntegrator::IntegrateSofteningCurve(double rArea, double rForceThreshold) const
{
    int iEnd = FindIndexJustAboveForceThreshold(rForceThreshold);
    double integral = IntegrateSofteningCurveInternal(iEnd);
    integral += IntegrateIndexToForceThreshold(iEnd, rForceThreshold);
    return integral / rArea;
}

double NuTo::Tools::GlobalFractureEnergyIntegrator::IntegrateIndexToForceThreshold(int rIndexEnd,
                                                                                   double rForceThreshold) const
{
    double d0 = mDispl(rIndexEnd - 1);
    double d1 = mDispl(rIndexEnd);
    double f0 = mForce(rIndexEnd - 1);
    double f1 = mForce(rIndexEnd);

    double dX = d0 - d1;
    double dF = f0 - f1;
    double dThreshold = d0 + dX / dF * (rForceThreshold - f0);
    // ToDo:: This interpolation should be done with NuTo::Math::Interpolation but I honestly have no idea how that
    // class works. TT.

    double meanForceStop = (rForceThreshold + f0) / 2.;
    double deltaDisplStop = dThreshold - d0;

    return meanForceStop * deltaDisplStop;
}

int NuTo::Tools::GlobalFractureEnergyIntegrator::FindIndexJustAboveForceThreshold(double rForceThreshold) const
{
    auto numRows = mForce.rows();
    if (rForceThreshold < mForce[numRows - 1])
        throw NuTo::Exception(__PRETTY_FUNCTION__, "The forces did not drop below the force threshold");

    for (auto i = 0; i < numRows - 1; ++i)
    {
        if (mForce[i] > rForceThreshold && mForce[i + 1] <= rForceThreshold)
            return i;
    }
    throw NuTo::Exception(__PRETTY_FUNCTION__, "The forces never reach the force threshold");
}

double NuTo::Tools::GlobalFractureEnergyIntegrator::IntegrateSofteningCurveInternal(int rEnd) const
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

void NuTo::Tools::GlobalFractureEnergyIntegrator::CheckData() const
{
    if (mForce.rows() != mDispl.rows())
        throw Exception(__PRETTY_FUNCTION__, "Force.Rows != Displacement.Rows()");

    if (std::abs(mForce[0]) > 1.e-10)
        throw Exception(__PRETTY_FUNCTION__, "Forces[0] != 0 ( but has to.)");

    if (std::abs(mDispl[0]) > 1.e-10)
        throw Exception(__PRETTY_FUNCTION__, "Displacements[0] != 0 ( but has to.)");
}
