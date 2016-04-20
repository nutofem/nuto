#include "nuto/math/BSplineCurve.h"


NuTo::BSplineCurve::BSplineCurve(int rDegree,
                                 const FullVector<double, Eigen::Dynamic>& rKnots,
                                 const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rControlPoints)
{
    int numKnots = rKnots.GetNumRows();
    int numControlPoints = rControlPoints.GetNumColumns();
    if (numControlPoints != numKnots - rDegree - 1)
            throw MathException("[BSplineCurve] incorrect number of control points or number of knots or degree.");

    mDegree = rDegree;
    mControlPoints = rControlPoints;
    mKnots = rKnots;
}

NuTo::BSplineCurve::BSplineCurve(int rDegree,
                                 const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints)
{
    assert(rDegree > 0);
    int numPoints = rPoints.GetNumRows();
    assert(numPoints);

    int dim = rPoints.GetNumColumns();

    mControlPoints.Resize(numPoints, dim);

    FullVector<double, Eigen::Dynamic> rParameters(numPoints);
    ParametrizationChordLengthMethod(rPoints, rParameters);

    // the coefficient matrix
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> A(numPoints, numPoints);

    FullVector<double, Eigen::Dynamic> basisFunctions(mDegree + 1);
    for (int i = 0; i < numPoints; i++)
    {
        int spanIdx = FindSpan(rParameters[i]);
        BasisFuns(rParameters[i], spanIdx, basisFunctions);
        A.SetBlock(i, spanIdx-mDegree, basisFunctions);
    }

    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> AInv (A.Inverse());

    mControlPoints = AInv*rPoints;


}

int NuTo::BSplineCurve::FindSpan(double rParameter)
{
    int numCP = GetNumControlPoints();
    if(rParameter == mKnots[numCP]) return numCP-1;

    int low  = mDegree;
    int high = numCP;
    int mid  = (low + high)/2;

    while(rParameter < mKnots[mid] || rParameter >= mKnots[mid+1])
    {
        if(rParameter < mKnots[mid]) high = mid;
        else                         low = mid;

        mid = (low + high)/2;
    }
    return mid;
}

void NuTo::BSplineCurve::BasisFuns(double rParameter, int rSpan, NuTo::FullVector<double, Eigen::Dynamic>& rBasisFunctions)
{
    assert(rBasisFunctions.GetNumRows() == (mDegree + 1));

    rBasisFunctions[0] = 1.;

    FullVector<double, Eigen::Dynamic> left(mDegree + 1);
    FullVector<double, Eigen::Dynamic> right(mDegree + 1);

    for (int j = 1; j <= mDegree; j++)
    {
        left[j]  = rParameter - mKnots[rSpan + 1 - j];
        right[j] = mKnots[rSpan + j] - rParameter;
        double saved = 0.;
        for(int r = 0; r < j; r++)
        {
            double temp = rBasisFunctions[r]/(right[r+1] + left[j-r]);
            rBasisFunctions[r] = saved + right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        rBasisFunctions[j] = saved;
    }
}

void NuTo::BSplineCurve::ParametrizationChordLengthMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                                                          FullVector<double, Eigen::Dynamic>& rParameters)
{
    int numPoints = rPoints.GetNumRows();
    assert(numPoints);
    rParameters.Resize(numPoints);

    rParameters[0] = 0.;
    rParameters[numPoints - 1] = 1.;

    double totalLengthPolygon = 0.;
    for (int i = 1; i <= numPoints - 1; i++)
    {
        NuTo::FullMatrix<double, 1, Eigen::Dynamic> diff(rPoints.GetRow(i) - rPoints.GetRow(i+1));
        totalLengthPolygon += diff.Norm();
    }

    for (int i = 1; i <= numPoints - 1; i++ )
    {
        NuTo::FullMatrix<double, 1, Eigen::Dynamic> diff(rPoints.GetRow(i) - rPoints.GetRow(i+1));
        rParameters[i] = rParameters[i-1] + diff.Norm()/totalLengthPolygon;
    }

    int numKnots = numPoints + mDegree;
    mKnots.resize(numKnots);

    for (int i = 0; i <= mDegree; i++) mKnots[i] = 0.;
    for (int i = numKnots - 1 - mDegree; i < numKnots; i++) mKnots[i] = 1.;

    double scale = 1./mDegree;
    for (int j = 1; j < (numKnots - 1 - mDegree); j++)
    {
        mKnots[j + mDegree] = 0.;
        for (int i = j; i <= (j + mDegree - 1); i++) mKnots[j + mDegree] += rParameters[i];
        mKnots[j + mDegree]/=scale;
    }
}


void NuTo::BSplineCurve::ParametrizationCentripetalMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints)
{

}
