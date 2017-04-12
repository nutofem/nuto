//
// Created by Thomas Titscher on 10/26/16.
//
#pragma once
#include <iostream>

#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"

namespace NuTo
{
namespace Test
{
    using namespace NuTo::Constitutive;
template <int TDim>
class ConstitutiveTangentTester
{
public:
    ConstitutiveTangentTester(IPConstitutiveLawBase& rLaw): mLaw(rLaw) {};
    ConstitutiveTangentTester(IPConstitutiveLawBase& rLaw, double rDelta, double rTolerance)
        : mLaw(rLaw), mDelta(rDelta), mTolerance(rTolerance) {}

    bool CheckTangent( ConstitutiveInputMap rInput, eInput rParameter, eOutput rValue, eOutput rTangent)
    {
        Eigen::MatrixXd tangent     = CalculateOutput(rInput, rTangent)->CopyToEigenMatrix();
        Eigen::MatrixXd tangent_cdf = CalculateTangentCDF(rInput, rParameter, rValue, rTangent);

        // If tangent is a ConstitutiveVector, its orientation is lost.
        // --> it is transposed, if needed.
        // Note that tangent_cdf HAS the right orientation
        if (tangent.rows() == tangent_cdf.cols() && tangent.cols() == tangent_cdf.rows())
            tangent.transposeInPlace();

        int row = 0, col = 0;
        Eigen::MatrixXd diffMatrix = (tangent - tangent_cdf).cwiseAbs();
        double maxEntry = diffMatrix.maxCoeff(&row, &col);

        std::cout << "Max. abs. error = " << maxEntry << " at (" << row << "," << col << ")   -  Tolerance = " << mTolerance << std::endl;
        if (maxEntry > mTolerance)
        {

            std::cout << "Tangent - Tangent_cdf: \n";
            std::cout << diffMatrix << "\n\n";

            std::cout << "Tangent: \n";
            std::cout << tangent << "\n\n";

            std::cout << "Tangent_cdf: \n";
            std::cout << tangent_cdf << "\n\n";
            return false;
        }

        return true;
    }

private:
    Eigen::MatrixXd CalculateTangentCDF(ConstitutiveInputMap rInput, eInput rParameter, eOutput rValue, eOutput rTangent)
    {
        auto value = CalculateOutput(rInput, rValue);
        auto& input = *rInput[rParameter];

        Eigen::MatrixXd tangent_cdf = Eigen::MatrixXd::Zero(value->GetNumRows(), input.GetNumRows());
        for (int i = 0; i < tangent_cdf.cols(); ++i)
        {
            input[i] -= mDelta / 2.;
            auto valueMinus = CalculateOutput(rInput, rValue);

            input[i] += mDelta;
            auto valuePlus = CalculateOutput(rInput, rValue);

            input[i] -= mDelta / 2; // back to normal

            for (int j = 0; j < tangent_cdf.rows(); ++j)
                tangent_cdf(i, j) = ( (*valuePlus)[j] - (*valueMinus)[j] ) / mDelta;
        }
        return tangent_cdf;
    }

    auto CalculateOutput(const ConstitutiveInputMap& rInput, eOutput rOutput)
    {
        ConstitutiveOutputMap output;
        output.Add<TDim>(rOutput);
        mLaw.Evaluate<TDim>(rInput, output);
        return output[rOutput]->clone();
    }

    IPConstitutiveLawBase& mLaw;
    double mDelta       = 1.e-6;
    double mTolerance   = 1.e-6;
};
}   // namespace Test
}   // namespace NuTo
