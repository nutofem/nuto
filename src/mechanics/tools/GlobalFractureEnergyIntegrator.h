//
// Created by Thomas Titscher on 10/28/16.
//
#pragma once
#include <string>
#include <Eigen/Core>

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
    GlobalFractureEnergyIntegrator(const Eigen::VectorXd& rForces, const Eigen::VectorXd& rDisplacements);

    //! @brief constructs from data files
    //! @param rForces <path>/<file>.dat to container with force data
    //! @param rDisplacements <path>/<file>.dat to container with force data
    //! @param rColumn column to read from, hint for most cases: 1=x, 2=y, 3=z direction
    //! @return obj
    GlobalFractureEnergyIntegrator(std::string rFileForces, std::string rFileDisplacements, int rColumn);

    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rArea "crack" area
    //! @return global GF
    double IntegrateSofteningCurve(double rArea) const;

    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rArea "crack" area
    //! @param rForceThreshold integration stops if the force drops below rForceThreshold
    //! @return global GF
    double IntegrateSofteningCurve(double rArea, double rForceThreshold) const;

private:
    //! @brief integrates F[rIndexEnd]
    //! @param rIndexEnd
    //! @param rForceThreshold
    //! @return
    double IntegrateIndexToForceThreshold(int rIndexEnd, double rForceThreshold) const;

    //! @brief finds the index where the force is just before the rForceThreshold
    //! @param rForceThreshold force threshold
    //! @return index of F[index] > rForceThreshold && F[index+1] < rForceThreshold
    int FindIndexJustAboveForceThreshold(double rForceThreshold) const;

    //! @brief calculates the global fracture energy by integrating over the load-displacement curve
    //! @param rEnd integrate from index=0 to index=rEnd
    //! @return global GF
    double IntegrateSofteningCurveInternal(int rEnd) const;

    void CheckData() const;

    Eigen::VectorXd mForce;
    Eigen::VectorXd mDispl;
};

} // namespace Tools
} // namespace NuTo
