#include "nuto/math/BSplineCurve.h"


NuTo::BSplineCurve::BSplineCurve(int rDegree,
                                 const FullVector<double, Eigen::Dynamic>& rKnots,
                                 const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rControlPoints)
{
    int numKnots = rKnots.GetNumRows();
    int numControlPoints = rControlPoints.GetNumColumns();

    if (numControlPoints != numKnots - rDegree - 1)
            throw MathException("[BSplineCurve] incorrect number of control points or number of knots or degree.");

    for (int i = 1; i < numKnots; i++)
        if(rKnots[i-1] > rKnots[i])
            throw MathException("[BSplineCurve] knot vector mus contain ascending values.");

    mDegree = rDegree;
    mControlPoints = rControlPoints;
    mKnots = rKnots;
}

NuTo::BSplineCurve::BSplineCurve(int rDegree,
                                 const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                                 FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& AInv)
{
    assert(rDegree > 0);
    mDegree = rDegree;
    int numPoints = rPoints.GetNumRows();
    assert(numPoints);

    int dim = rPoints.GetNumColumns();

    mControlPoints.Resize(numPoints, dim);

    FullVector<double, Eigen::Dynamic> rParameters = ParametrizationChordLengthMethod(rPoints);
    ParametrizationKnotVector(rParameters);

    // the coefficient matrix
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> A(numPoints, numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        int spanIdx = FindSpan(rParameters[i]);
        FullMatrix<double, 1, Eigen::Dynamic> basisFunctions = BasisFuns(rParameters[i], spanIdx);
        A.SetBlock(i, spanIdx-mDegree, basisFunctions);
    }

    AInv = A.Inverse();

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

NuTo::FullMatrix<double, 1, Eigen::Dynamic> NuTo::BSplineCurve::BasisFuns(double rParameter, int rSpan)
{
    FullMatrix<double, 1, Eigen::Dynamic> rBasisFunctions(1, mDegree + 1);

    rBasisFunctions.SetValue(0,0,1.);

    FullVector<double, Eigen::Dynamic> left(mDegree + 1);
    FullVector<double, Eigen::Dynamic> right(mDegree + 1);

    for (int j = 1; j <= mDegree; j++)
    {
        left[j]  = rParameter - mKnots[rSpan + 1 - j];
        right[j] = mKnots[rSpan + j] - rParameter;
        double saved = 0.;
        for(int r = 0; r < j; r++)
        {
            double temp = rBasisFunctions(0,r)/(right[r+1] + left[j-r]);
            rBasisFunctions(0,r) = saved + right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        rBasisFunctions(0,j) = saved;
    }

    return rBasisFunctions;
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::BSplineCurve::BasisFunsAndDerivatives(double rParameter, int rSpan, int maxDer)
{
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ndu(mDegree, mDegree);

    ndu(0,0) = 1.0;

    FullVector<double, Eigen::Dynamic> left(mDegree + 1);
    FullVector<double, Eigen::Dynamic> right(mDegree + 1);

    for (int j = 1; j <= mDegree; j++)
    {
        left[j]  = rParameter - mKnots[rSpan + 1 - j];
        right[j] = mKnots[rSpan + j] - rParameter;
        double saved = 0.;
        for(int r = 0; r < j; r++)
        {
            ndu(j,r) = right[r+1] + left[j-r]; // lower triange (knot differences)
            double temp = ndu(r, j-1)/ndu(j,r);
            ndu(r,j) = saved + right[r+1]*temp; // upper triangle
            saved = left[j-r]*temp;
        }
        ndu(j,j) = saved;
    }

    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ders(mDegree, mDegree);

    for(int j = 0; j <= mDegree; j++)
    {
        ders(0,j) = ndu(j,mDegree);

        FullMatrix<double, 2, 2> a(2, mDegree+1);
        int r = 0;
        for(; r<=mDegree; r++)
        {
            int s1 = 0;
            int s2 = 1;

            a(0,0) = 1.0;

            for(int k = 1; k<=maxDer; k++)
            {
                double d = 0.0;
                int rk = r-k;
                int pk = mDegree-k;
                if(r >= k)
                {
                    a(s2,0) = a(s1,0)/ndu(pk+1, rk);
                    double d = a(s2,0)*ndu(rk,pk);
                }
                int j1 = 0;
                if(rk >= -1) j1 = 1;
                else         j1 = -rk;

                int j2 = 0;
                if(r-1 <= pk) j2 = k-1;
                else          j2 = mDegree-r;

                for(j = j1; j <= j2; j++)
                {
                    a(s2, j) = (a(s1,j) - a(s1, j-1))/ndu(pk+1, rk+j);
                    d += a(s2,j)*ndu(rk+j, pk);
                }

                if(r <= pk)
                {
                    a(s2,k) = -a(s1,k-1)/ndu(pk+1,r);
                    d += a(s2,k)*ndu(r,pk);
                }

                ders(k,r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
            }
        }

        r = mDegree;
        for(int k = 1; k <= maxDer; k++)
        {
            for (int j = 0; j <= mDegree; j++) ders(k,j) *=r;
            r*= (mDegree-k);
        }
    }

    return ders.Trans();

}


NuTo::FullVector<double, Eigen::Dynamic> NuTo::BSplineCurve::ParametrizationCentripetalMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints)
{
    int numPoints = rPoints.GetNumRows();
    assert(numPoints);

    FullVector<double, Eigen::Dynamic> rParameters;
    rParameters.Resize(numPoints);

    rParameters[0] = 0.;
    rParameters[numPoints - 1] = 1.;

    double totalLengthPolygon = 0.;
    for (int i = 1; i <= (numPoints-1); i++)
    {
        NuTo::FullMatrix<double, 1, Eigen::Dynamic> diff(rPoints.GetRow(i) - rPoints.GetRow(i-1));
        totalLengthPolygon += std::sqrt(diff.Norm());
    }

    for (int i = 1; i < (numPoints-1); i++ )
    {
        NuTo::FullMatrix<double, 1, Eigen::Dynamic> diff(rPoints.GetRow(i) - rPoints.GetRow(i-1));
        rParameters[i] = rParameters[i-1] + std::sqrt(diff.Norm())/totalLengthPolygon;
    }

    return rParameters;
}

NuTo::FullVector<double, Eigen::Dynamic> NuTo::BSplineCurve::ParametrizationChordLengthMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints)
{
    int numPoints = rPoints.GetNumRows();
    assert(numPoints);

    FullVector<double, Eigen::Dynamic> rParameters;
    rParameters.Resize(numPoints);

    rParameters[0] = 0.;
    rParameters[numPoints - 1] = 1.;

    double totalLengthPolygon = 0.;
    for (int i = 1; i <= (numPoints-1); i++)
    {
        NuTo::FullMatrix<double, 1, Eigen::Dynamic> diff(rPoints.GetRow(i) - rPoints.GetRow(i-1));
        totalLengthPolygon += diff.Norm();
    }

    for (int i = 1; i < (numPoints-1); i++ )
    {
        NuTo::FullMatrix<double, 1, Eigen::Dynamic> diff(rPoints.GetRow(i) - rPoints.GetRow(i-1));
        rParameters[i] = rParameters[i-1] + diff.Norm()/totalLengthPolygon;
    }

    return rParameters;
}

const NuTo::FullVector<double, Eigen::Dynamic>& NuTo::BSplineCurve::ParametrizationKnotVector(const FullVector<double, Eigen::Dynamic>& rParameters)
{
    int numPoints = rParameters.GetNumRows();
    assert(numPoints);
    int numKnots = numPoints + mDegree + 1;
    mKnots.resize(numKnots);

    for (int i = 0; i <= mDegree; i++) mKnots[i] = 0.;
    for (int i = numKnots - 1 - mDegree; i < numKnots; i++) mKnots[i] = 1.;

    double scale = 1./mDegree;
    for (int j = 1; j <= (numPoints - 1 - mDegree); j++)
    {
        mKnots[j + mDegree] = 0.;
        for (int i = j; i <= (j + mDegree - 1); i++) mKnots[j + mDegree] += rParameters[i];
        mKnots[j + mDegree]*=scale;
    }

    return mKnots;
}

NuTo::FullVector<double, Eigen::Dynamic> NuTo::BSplineCurve::CurvePoint(double rParameter)
{
    int spanIdx = FindSpan(rParameter);
    FullMatrix<double, 1, Eigen::Dynamic> basisFunctions = BasisFuns(rParameter, spanIdx);

    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates(1, GetDimension());
    for(int i = 0; i <= mDegree; i++) coordinates += basisFunctions.GetValue(0,i)*(mControlPoints.GetRow(spanIdx-mDegree+i));

    return coordinates.Trans();
}


NuTo::FullMatrix<double, Eigen::Dynamic> NuTo::BSplineCurve::CurvePoints(FullVector<double, Eigen::Dynamic> rParameter)
{
    int numParameters = rParameter.GetNumRows();
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> result(numParameters, GetDimension());
    for(int i = 0; i < numParameters; i++)
    {
        result.SetBlock(i,0,CurvePoint(rParameter[i]).Trans());
    }
    return result;
}
