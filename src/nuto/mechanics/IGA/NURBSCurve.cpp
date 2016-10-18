#include "nuto/mechanics/IGA/NURBSCurve.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"

NuTo::NURBSCurve::NURBSCurve(const Eigen::MatrixXd &rKnots,
                                 const Eigen::MatrixXd &rControlPoints,
                                 const Eigen::VectorXd &rWeights,
                                 int rDegree)
{
    int numKnots = rKnots.rows();
    int numControlPoints = rControlPoints.rows();

    if (numControlPoints != numKnots - 1 - rDegree)
        throw Exception("[NURBSCurve] incorrect number of control points or number of knots or degree.");

    for (int i = 1; i < numKnots; i++)
        if(rKnots(i-1) > rKnots(i))
            throw Exception("[NURBSCurve] knot vector must contain ascending values.");

    mDegree = rDegree;
    mControlPoints = rControlPoints;
    mKnots = rKnots;
    mWeights = rWeights;
}

NuTo::NURBSCurve::NURBSCurve(int rDegree,
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

    mWeights.resize(mControlPoints.rows(), mControlPoints.cols());
    mWeights.setOnes(mControlPoints.rows(), mControlPoints.cols());

    // the coefficient matrix
    Eigen::MatrixXd A(numPoints, numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        int spanIdx = ShapeFunctionsIGA::FindSpan(rParameters[i], mDegree, mKnots);
        Eigen::VectorXd basisFunctions = ShapeFunctionsIGA::BasisFunctionsRat(rParameters[i], spanIdx, mDegree, mKnots, mWeights);
        A.block(i, spanIdx-mDegree, 1, basisFunctions.rows()) = basisFunctions.transpose();
    }


    AInv = A.inverse();

    mControlPoints = AInv*rPoints;
}

int NuTo::NURBSCurve::GetNumIGAElements() const
{
    int numElements = 0;

    for(int i = 1; i < mKnots.rows(); i++)
    {
        if(mKnots(i-1) < mKnots(i)) numElements++;
    }

    return numElements;
}

int NuTo::NURBSCurve::GetElementFirstKnotID(int rElementID) const
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
        throw Exception("[NURBSCurve::GetElementFirstKnotID] ElementId not found.");

    return knotID;
}

Eigen::MatrixXd NuTo::NURBSCurve::GetElementKnots(int rElementID) const
{
    int knotID = GetElementFirstKnotID(rElementID);
    Eigen::MatrixXd knots(1,2);

    knots << mKnots(knotID) , mKnots(knotID+1);

    return knots;
}

Eigen::VectorXi NuTo::NURBSCurve::GetElementKnotIDs(int rElementID) const
{
    int knotID = GetElementFirstKnotID(rElementID);
    Eigen::VectorXi knots(1);

    knots << knotID ;

    return knots;
}

Eigen::MatrixXd NuTo::NURBSCurve::GetElementControlPoints(int rElementID) const
{
    assert((rElementID >= 0) && (rElementID < GetNumIGAElements()));

    int knotID = GetElementFirstKnotID(rElementID);

    int rows = mControlPoints.rows();

    return mControlPoints.block(0, knotID - mDegree, rows, mDegree+1);
}

Eigen::VectorXi NuTo::NURBSCurve::GetElementControlPointIDs(int rElementID) const
{
    assert((rElementID >= 0) && (rElementID < GetNumIGAElements()));

    int knotID = GetElementFirstKnotID(rElementID);

    Eigen::VectorXi ids(mDegree+1);

    for(int i = 0; i < mDegree + 1; i++) ids(i) = knotID - mDegree + i;

    return ids;
}

int NuTo::NURBSCurve::GetMultiplicityOfKnot(double rKnot)
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rKnot, mDegree, mKnots);

    int count = 0;
    for(int i = spanIdx; i > 0; i--)
        if(mKnots(i) == rKnot) count++;
        else                             break;

    return count;
}

const Eigen::MatrixXd& NuTo::NURBSCurve::GetBezierExtraction(int rElementID) const
{
    assert(((size_t)rElementID > (mBezierOperators.size()-1)) || (rElementID < 0));
    if (mBezierOperators.size() == 0)
        throw Exception("[NURBSCurve] The Bézier extraction has to be calculated first. Use the function 'BezierElementExtractionOperators'.");

    return mBezierOperators[rElementID];
}

const std::vector<Eigen::MatrixXd>& NuTo::NURBSCurve::GetBezierExtraction() const
{
    if (mBezierOperators.size() == 0)
        throw Exception("[NURBSCurve] The Bézier extraction has to be calculated first. Use the function 'BezierElementExtractionOperators'.");

    return mBezierOperators;
}

Eigen::VectorXd NuTo::NURBSCurve::BasisFunctions(double rParameter) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    return ShapeFunctionsIGA::BasisFunctionsRat(rParameter, spanIdx, mDegree, mKnots, mWeights);
}

Eigen::MatrixXd NuTo::NURBSCurve::BasisFunctionsAndDerivatives(double rParameter, int der) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    return ShapeFunctionsIGA::BasisFunctionsAndDerivativesRat(der, rParameter, spanIdx, mDegree, mKnots, mWeights);
}

Eigen::VectorXd NuTo::NURBSCurve::ParametrizationCentripetalMethod(const Eigen::MatrixXd& rPoints) const
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

Eigen::VectorXd NuTo::NURBSCurve::ParametrizationChordLengthMethod(const Eigen::MatrixXd& rPoints) const
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

const Eigen::VectorXd& NuTo::NURBSCurve::ParametrizationKnotVector(const Eigen::VectorXd& rParameters)
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

void NuTo::NURBSCurve::BezierElementExtractionOperators()
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

void NuTo::NURBSCurve::BezierElementControlPoints()
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

Eigen::VectorXd NuTo::NURBSCurve::CurvePoint(double rParameter) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    Eigen::VectorXd basisFunctions = ShapeFunctionsIGA::BasisFunctionsRat(rParameter, spanIdx, mDegree, mKnots, mWeights);

    Eigen::VectorXd coordinates(GetDimension());
    coordinates.setZero(GetDimension());

    for(int i = 0; i < GetDimension(); i++)
        for(int j = 0; j <= mDegree; j++)
            coordinates(i) += basisFunctions(j)*mControlPoints(spanIdx-mDegree+j, i);

    return coordinates;
}


Eigen::MatrixXd NuTo::NURBSCurve::CurvePoints(const Eigen::VectorXd &rParameter) const
{
    int numParameters = rParameter.rows();
    Eigen::MatrixXd result(numParameters, GetDimension());
    for(int i = 0; i < numParameters; i++)
    {
        result.row(i) = CurvePoint(rParameter[i]).transpose();
    }

    return result;
}

void NuTo::NURBSCurve::InsertKnot(double rKnotToInsert, int rMultiplicity)
{
    int k = ShapeFunctionsIGA::FindSpan(rKnotToInsert, mDegree, mKnots);
    int initialMultiplicity = GetMultiplicityOfKnot(rKnotToInsert);
    assert(initialMultiplicity + rMultiplicity <= mDegree);
    // new knot vector
    Eigen::VectorXd newKnots(GetNumKnots() + rMultiplicity);

    for(int i = 0 ; i <= k; i++)               newKnots(i) = mKnots(i);
    for(int i = 1 ; i <= rMultiplicity; i++)   newKnots(k + i) = rKnotToInsert;
    for(int i = k + 1; i < GetNumKnots(); i++) newKnots(rMultiplicity + i) = mKnots(i);

    // control points
    Eigen::MatrixXd newControlPoints(GetNumControlPoints() + rMultiplicity, GetDimension());
    Eigen::MatrixXd Rw(mDegree+1, GetDimension());

    for(int i = 0; i <= k - mDegree; i++) newControlPoints.row(i) = mControlPoints.row(i);
    for(int i = k - initialMultiplicity; i < GetNumControlPoints(); i++) newControlPoints.row(i + rMultiplicity) = mControlPoints.row(i);
    for(int i = 0; i <= mDegree - initialMultiplicity; i++) Rw.row(i) = mControlPoints.row(k - mDegree + i);

    int L = 0;
    for(int j = 1; j <= rMultiplicity; j++)
    {
        L = k - mDegree + j;
        for(int i = 0; i <= mDegree - j -initialMultiplicity; i++)
        {
            double alpha = (rKnotToInsert - mKnots(L+i))/(mKnots(i+k+1) - mKnots(L+i));
            Rw.row(i) = alpha*Rw.row(i+1) + (1. - alpha)*Rw.row(i);
        }
        newControlPoints.row(L) = Rw.row(0);
        newControlPoints.row(k + rMultiplicity - j - initialMultiplicity) = Rw.row(mDegree - j - initialMultiplicity);
    }

    for(int i = L + 1; i < k - initialMultiplicity; i++) newControlPoints.row(i) = Rw.row(i-L);

    mKnots = newKnots;
    mControlPoints = newControlPoints;
}


void NuTo::NURBSCurve::RefineKnots(const Eigen::VectorXd &rKnotsToInsert)
{
    int numInsert = rKnotsToInsert.rows();
    int begin     = ShapeFunctionsIGA::FindSpan(rKnotsToInsert(0),           mDegree, mKnots);
    int end       = ShapeFunctionsIGA::FindSpan(rKnotsToInsert(numInsert-1), mDegree, mKnots);
    end++;

    // since we're dealing with NURBS curves, a projection is needed ...
    Eigen::MatrixXd controlPointsProjected(mControlPoints.rows(), GetDimension()+1);
    controlPointsProjected.block(0, 0, mControlPoints.rows(), mControlPoints.cols()) = mControlPoints;

    int col = mControlPoints.cols();
    for(int i = 0; i <  mControlPoints.rows(); i++)
        controlPointsProjected(i, col) = mWeights(i);

    // new control points
    Eigen::MatrixXd newControlPoints(GetNumControlPoints() + numInsert, GetDimension()+1);
    for(int i = 0; i <= begin - mDegree; i++)  newControlPoints.row(i) = controlPointsProjected.row(i);
    for(int i = end - 1; i < GetNumControlPoints(); i++) newControlPoints.row(i + numInsert) = controlPointsProjected.row(i);

    // new knot vector
    Eigen::VectorXd newKnots(GetNumKnots() + numInsert);
    for(int i = 0 ; i <= begin; i++) newKnots(i) = mKnots(i);
    for(int i = end + mDegree; i < GetNumKnots(); i++) newKnots(numInsert + i) = mKnots(i);

    int i = end + mDegree - 1;
    int k = end + mDegree + numInsert - 1;

    for(int j = numInsert - 1; j >= 0; j--)
    {
        while(rKnotsToInsert(j) <= mKnots(i) && i > begin)
        {
            newControlPoints.row(k - mDegree - 1) = controlPointsProjected.row(i - mDegree - 1);
            newKnots(k) = mKnots(i);
            k--;
            i--;
        }
        newControlPoints.row(k - mDegree - 1) = newControlPoints.row(k - mDegree);
        for(int l = 1; l <= mDegree; l++)
        {
            int ind = k - mDegree + l;
            double alpha = newKnots(k + l) - rKnotsToInsert(j);
            if(std::fabs(alpha) == 0.0)
            {
                newControlPoints.row(ind -1) = newControlPoints.row(ind);
            }
            else
            {
                alpha = alpha/(newKnots(k+l) - mKnots(i - mDegree + l));
                newControlPoints.row(ind - 1) = alpha*newControlPoints.row(ind - 1) + (1.0 - alpha)*newControlPoints.row(ind);
            }
        }
        newKnots(k) = rKnotsToInsert(j);
        k--;
    }

    mKnots = newKnots;
    mWeights = newControlPoints.col(newControlPoints.cols() - 1);
    mControlPoints.resize(newControlPoints.rows(), GetDimension());
    for(int i = 0; i < newControlPoints.rows(); i++) mControlPoints.row(i) = newControlPoints.block(i,0, 1, GetDimension())/mWeights(i);
}

void NuTo::NURBSCurve::DuplicateKnots()
{
    Eigen::VectorXd knotsToInsert(GetNumIGAElements());
    for (int i = 0; i < knotsToInsert.rows(); i++)
    {
        int beginID = GetElementFirstKnotID(i);
        knotsToInsert(i) = (mKnots(beginID + 1) + mKnots(beginID))/2.;
    }
    RefineKnots(knotsToInsert);
}
