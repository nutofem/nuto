#include "nuto/mechanics/IGA/BSplineCurve.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"

NuTo::BSplineCurve::BSplineCurve(const Eigen::MatrixXd& rKnots,
                                 const Eigen::MatrixXd& rControlPoints,
                                 int rDegree)
{
    int numKnots = rKnots.rows();
    int numControlPoints = rControlPoints.rows();

    if (numControlPoints != numKnots - 1 - rDegree)
            throw Exception("[BSplineCurve] incorrect number of control points or number of knots or degree.");

    for (int i = 1; i < numKnots; i++)
        if(rKnots(i-1) > rKnots(i))
            throw Exception("[BSplineCurve] knot vector must contain ascending values.");

    mDegree = rDegree;
    mControlPoints = rControlPoints;
    mKnots = rKnots;
}

NuTo::BSplineCurve::BSplineCurve(int rDegree,
                                 const Eigen::MatrixXd& rPoints,
                                 Eigen::MatrixXd& AInv)
{
    assert(rDegree > 0);
    mDegree = rDegree;
    int numPoints = rPoints.rows();
    assert(numPoints);

    int dim = rPoints.cols();

    mControlPoints.resize(numPoints, dim);

    Eigen::VectorXd rParameters = ParametrizationChordLengthMethod(rPoints);
    ParametrizationKnotVector(rParameters);

    // the coefficient matrix
    Eigen::MatrixXd A(numPoints, numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        int spanIdx = ShapeFunctionsIGA::FindSpan(rParameters[i], mDegree, mKnots);
        Eigen::VectorXd basisFunctions = ShapeFunctionsIGA::BasisFunctions(rParameters[i], spanIdx, mDegree, mKnots);
        A.block(i, spanIdx-mDegree, 1, basisFunctions.rows()) = basisFunctions.transpose();
    }

    AInv = A.inverse();

    mControlPoints = AInv*rPoints;

}

int NuTo::BSplineCurve::GetNumIGAElements() const
{
    int numElements = 0;

    for(int i = 1; i < mKnots.rows(); i++)
    {
        if(mKnots(i-1) < mKnots(i)) numElements++;
    }

    return numElements;
}

int NuTo::BSplineCurve::GetElementFirstKnotID(int rElementID) const
{
    int elementID = 0;
    int knotID = -1;
    for(int i = 1; i < mKnots.rows(); i++)
    {
        if(mKnots(i-1) < mKnots(i)) elementID++;
        if(rElementID == elementID-1)
        {
            knotID = i-1;
            break;
        }
    }
    if (knotID == -1)
        throw Exception("[BSplineCurve::GetElementFirstKnotID] ElementId not found.");

    return knotID;
}

Eigen::MatrixXd NuTo::BSplineCurve::GetElementKnots(int rElementID) const
{
    int knotID = GetElementFirstKnotID(rElementID);
    Eigen::MatrixXd knots(1,2);

    knots << mKnots(knotID) , mKnots(knotID+1);

    return knots;
}

Eigen::VectorXi NuTo::BSplineCurve::GetElementKnotIDs(int rElementID) const
{
    int knotID = GetElementFirstKnotID(rElementID);
    Eigen::VectorXi knots(1);

    knots << knotID ;

    return knots;
}

Eigen::MatrixXd NuTo::BSplineCurve::GetElementControlPoints(int rElementID) const
{
    int numIGAElements = GetNumIGAElements();
    assert((rElementID >= 0) && (rElementID < numIGAElements));

    int knotID = GetElementFirstKnotID(rElementID);

    int rows = mControlPoints.rows();

    return mControlPoints.block(0, knotID - mDegree, rows, mDegree+1);
}

Eigen::VectorXi NuTo::BSplineCurve::GetElementControlPointIDs(int rElementID) const
{
    int numIGAElements = GetNumIGAElements();
    assert((rElementID >= 0) && (rElementID < numIGAElements));

    int knotID = GetElementFirstKnotID(rElementID);

    Eigen::VectorXi ids(mDegree+1);

    for(int i = 0; i < mDegree + 1; i++) ids(i) = knotID - mDegree + i;

    return ids;
}

const Eigen::MatrixXd& NuTo::BSplineCurve::GetBezierExtraction(int rElementID) const
{
    int size = mBezierOperators.size();
    assert((rElementID > (size-1)) || (rElementID < 0));
    if (mBezierOperators.size() == 0)
        throw Exception("[BSplineCurve] The Bézier extraction has to be calculated first. Use the function 'BezierElementExtractionOperators'.");

    return mBezierOperators[rElementID];
}

const std::vector<Eigen::MatrixXd>& NuTo::BSplineCurve::GetBezierExtraction() const
{
    if (mBezierOperators.size() == 0)
        throw Exception("[BSplineCurve] The Bézier extraction has to be calculated first. Use the function 'BezierElementExtractionOperators'.");

    return mBezierOperators;
}

Eigen::VectorXd NuTo::BSplineCurve::BasisFunctions(double rParameter) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    return ShapeFunctionsIGA::BasisFunctions(rParameter, spanIdx, mDegree, mKnots);
}

Eigen::MatrixXd NuTo::BSplineCurve::BasisFunctionsAndDerivatives(double rParameter, int maxDer) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    return ShapeFunctionsIGA::BasisFunctionsAndDerivatives(rParameter, spanIdx, maxDer, mDegree, mKnots);
}

Eigen::VectorXd NuTo::BSplineCurve::ParametrizationCentripetalMethod(const Eigen::MatrixXd& rPoints) const
{
    int numPoints = rPoints.rows();
    assert(numPoints);

    Eigen::VectorXd rParameters;
    rParameters.resize(numPoints);

    rParameters[0] = 0.;
    rParameters[numPoints - 1] = 1.;

    double totalLengthPolygon = 0.;
    for (int i = 1; i <= (numPoints-1); i++)
    {
        Eigen::Matrix<double, 1, Eigen::Dynamic> diff(rPoints.row(i) - rPoints.row(i-1));
        totalLengthPolygon += std::sqrt(diff.norm());
    }

    for (int i = 1; i < (numPoints-1); i++ )
    {
        Eigen::Matrix<double, 1, Eigen::Dynamic> diff(rPoints.row(i) - rPoints.row(i-1));
        rParameters[i] = rParameters[i-1] + std::sqrt(diff.norm())/totalLengthPolygon;
    }

    return rParameters;
}

Eigen::VectorXd NuTo::BSplineCurve::ParametrizationChordLengthMethod(const Eigen::MatrixXd& rPoints) const
{
    int numPoints = rPoints.rows();
    assert(numPoints);

    Eigen::VectorXd rParameters;
    rParameters.resize(numPoints);

    rParameters[0] = 0.;
    rParameters[numPoints - 1] = 1.;

    double totalLengthPolygon = 0.;
    for (int i = 1; i <= (numPoints-1); i++)
    {
        Eigen::MatrixXd diff(rPoints.row(i) - rPoints.row(i-1));
        totalLengthPolygon += diff.norm();
    }

    for (int i = 1; i < (numPoints-1); i++ )
    {
        Eigen::MatrixXd diff(rPoints.row(i) - rPoints.row(i-1));
        rParameters[i] = rParameters[i-1] + diff.norm()/totalLengthPolygon;
    }

    return rParameters;
}

const Eigen::VectorXd& NuTo::BSplineCurve::ParametrizationKnotVector(const Eigen::VectorXd& rParameters)
{
    int numPoints = rParameters.rows();
    assert(numPoints);
    int numKnots = numPoints + mDegree + 1;
    mKnots.resize(numKnots);

    for (int i = 0; i <= mDegree; i++) mKnots(i) = 0.;
    for (int i = numKnots - 1 - mDegree; i < numKnots; i++) mKnots(i) = 1.;

    double scale = 1./mDegree;
    for (int j = 1; j <= (numPoints - 1 - mDegree); j++)
    {
        mKnots(j + mDegree) = 0.;
        for (int i = j; i <= (j + mDegree - 1); i++) mKnots(j + mDegree) += rParameters[i];
        mKnots(j + mDegree)*=scale;
    }

    return mKnots;
}

void NuTo::BSplineCurve::BezierElementExtractionOperators()
{
    int formerKnotIndex = mDegree; // 0 ... mDegree - equal knots for endpoint interpolation
    int numElements = 1;
    int numKnots = GetNumKnots();

    mBezierOperators.resize(GetNumIGAElements());
    for(std::vector<Eigen::MatrixXd>::iterator it = mBezierOperators.begin(); it != mBezierOperators.end(); it++)
    {
        (*it).resize(mDegree+1, mDegree+1);
        (*it).setIdentity();
    }

    for(int currentKnotIndex = formerKnotIndex + 1; currentKnotIndex < numKnots; )
    {
        double knotToInsert = mKnots(currentKnotIndex);

        int i = currentKnotIndex;
        for(; (currentKnotIndex < numKnots-1) && (mKnots(currentKnotIndex+1) == mKnots(currentKnotIndex)); currentKnotIndex++);
        int mult = currentKnotIndex-i+1;

        if(mult < mDegree)
        {
            Eigen::VectorXd alphas(numKnots);
            for(int A = 0; A < currentKnotIndex - mDegree; A++) alphas(A) = 1.;
            for(int A = currentKnotIndex - mDegree ; A < currentKnotIndex; A++) alphas(A) = (knotToInsert - mKnots(A))/(mKnots(A+mDegree) - mKnots(A));
            for(int A = currentKnotIndex; A < numKnots; A++) alphas(A) = 0.;

            double diff = mKnots(currentKnotIndex) - mKnots(formerKnotIndex);
            for(int j = mDegree ; j <= mult+1; j++)
            {
                double temp =  diff/(mKnots(formerKnotIndex+j) - mKnots(formerKnotIndex));
                alphas(j-mult-1) = temp;
            }

            int r = mDegree - mult;

            for(int j = 1; j <= r; j++)
            {
                int save = r-j+1;
                int s = mult+j;

                for(int k = mDegree + 1; k <= s+1; k++)
                {
                    double alpha = alphas(k-s-1);
                    mBezierOperators[numElements-1].col(k-1) = alpha*mBezierOperators[numElements-1].col(k-1) + (1. - alpha)*mBezierOperators[numElements-1].col(k-2);
                }
                if(currentKnotIndex < numKnots)
                {
                    mBezierOperators[numElements].block(save-1, save-1, j+1, 1) = mBezierOperators[numElements-1].block(mDegree-j, mDegree, j+1, 1);
                }
            }
            numElements++;
        }
        formerKnotIndex=currentKnotIndex;
        currentKnotIndex++;
    }
}

void NuTo::BSplineCurve::BezierElementControlPoints()
{
    mBezierControlPoints.resize(GetNumIGAElements());

    int dim = mControlPoints.cols();

    std::vector<Eigen::MatrixXd>::iterator itCP = mBezierControlPoints.begin();
    int i = 0;
    for(std::vector<Eigen::MatrixXd>::iterator it = mBezierOperators.begin(); it != mBezierOperators.end(); it++, itCP++, i++)
    {
        Eigen::MatrixXd controlPoints = (mControlPoints.block(0, i, mDegree+1, dim)).transpose();
        (*itCP) = (*it)*controlPoints;
    }
}

Eigen::VectorXd NuTo::BSplineCurve::CurvePoint(double rParameter) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    Eigen::VectorXd basisFunctions = ShapeFunctionsIGA::BasisFunctions(rParameter, spanIdx, mDegree, mKnots);

    Eigen::VectorXd coordinates(GetDimension());
    for(int i = 0; i <= mDegree; i++) coordinates += basisFunctions(i)*(mControlPoints.row(spanIdx-mDegree+i));

    return coordinates;
}


Eigen::MatrixXd NuTo::BSplineCurve::CurvePoints(const Eigen::VectorXd &rParameter) const
{
    int numParameters = rParameter.rows();
    Eigen::MatrixXd result(GetDimension(), numParameters);
    for(int i = 0; i < numParameters; i++)
    {
        result.col(i) = CurvePoint(rParameter[i]);
    }
    return result;
}
