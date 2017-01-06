#include "math/EigenCompanion.h"
#include "math/MathException.h"

using namespace NuTo;

void EigenCompanion::AppendRows(Eigen::MatrixXd& top, const Eigen::MatrixXd& bottom)
{
	if (top.cols() != bottom.cols())
		throw MathException(__PRETTY_FUNCTION__, "Number of columns for both matrices must be identical.");
    top.conservativeResize(top.rows() + bottom.rows(), Eigen::NoChange);
    top.bottomRows(bottom.rows()) = bottom;
}
