#include "nuto/mechanics/tools/ImplicitTimeIntegration.h"

#include <ostream>
#include <iomanip>
#include <boost/range/numeric.hpp>

#include "nuto/base/Logger.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"

using namespace NuTo;

//// --------------------- Implicit time integration --------------------- //
// void ImplicitTimeIntegration::WriteTimeDofResidual(Eigen::VectorXd& u, std::ostream& out, QuasiStaticProblem&
// callBack)
//{
//    /* Each call to ToDofVector<double> assumes that we solve for the step t + dt. But we are in postprocess right
//     * now.
//     * This we are already updated to t + dt. Keeping mGlobalTime as it is, will cause the ToDofVector<double> to
//     * apply
//     * constraints for t + dt + dt.
//     * keep in mind that the gradient is not the residual, since constraint dofs do not have a non-vanishing gradient
//     */

//    auto residual = callBack.FullResidual(u);
//    double residualSum = 0; // boost::accumulate(residual, 0.);

//    out << callBack.GetGlobalTime() << '\t' << residualSum << '\n';
//    out << std::flush; // We really want to flush the output in case the program is interupted by CTRL-C or
//    // something.
//}


int ImplicitTimeIntegration::DoStep(DofVector<double>& start, QuasiStaticProblem& callBack, double newGlobalTime,
                                    std::string solverType)
{
    // compute trial solution (includes update of the constraint dofs, no line search)
    // trialU is with constrains
    Eigen::VectorXd trialU = callBack.TrialState(start, newGlobalTime, solverType);

    int numIterations = 0;

    Eigen::VectorXd tmpX;
    try
    {
        tmpX = NewtonRaphson::Solve(callBack, trialU, EigenSparseSolver(solverType), 6, NewtonRaphson::LineSearch(),
                                    &numIterations);
    }
    catch (std::exception& e)
    {
        throw NewtonRaphson::NoConvergence(e.what());
    }

    if (callBack.Norm(tmpX) > 1.e10)
        throw NewtonRaphson::NoConvergence("", "floating point exception");

    callBack.UpdateHistory(tmpX);
    callBack.ToDofVector2(tmpX, start);

    return numIterations;
}
