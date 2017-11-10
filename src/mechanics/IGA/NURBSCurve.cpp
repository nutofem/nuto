/*
 * NURBSCurve.cpp
 *
 *  Created on: Nov 2, 2017
 *      Author: vkindrac: taken from Peter Otto
 */

#include "mechanics/IGA/NURBSCurve.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

NuTo::NURBSCurve::NURBSCurve(const Eigen::VectorXd &rKnots,
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

NuTo::NURBSCurve::NURBSCurve(int rDegree, const Eigen::MatrixXd &rPoints, Eigen::MatrixXd &A, NuTo::NURBSCurve::mParametrization rParametrizationMethod)
{
    assert(rDegree > 0);
    mDegree = rDegree;
    int numPoints = rPoints.rows();
    assert(numPoints >= mDegree);

    int dim = rPoints.cols();

    mControlPoints.resize(numPoints, dim);

    Eigen::VectorXd rParameters;
    switch (rParametrizationMethod) {
    case mParametrization::chord:
        rParameters = ParametrizationChordLengthMethod(rPoints);
        break;
    case mParametrization::centripetal:
        rParameters = ParametrizationCentripetalMethod(rPoints);
        break;
    default:
        break;
    }

    ParametrizationKnotVector(rParameters);

    mWeights.resize(numPoints, 1);
    mWeights.setOnes(numPoints, 1);

    // the coefficient matrix
    A.resize(numPoints, numPoints);
    A.setZero(numPoints, numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        int spanIdx = ShapeFunctionsIGA::FindSpan(rParameters[i], mDegree, mKnots);
        Eigen::VectorXd basisFunctions = ShapeFunctionsIGA::BasisFunctionsRat(rParameters[i], spanIdx, mDegree, mKnots, mWeights);
        A.block(i, spanIdx-mDegree, 1, basisFunctions.rows()) = basisFunctions.transpose();
    }

    mControlPoints = A.fullPivLu().solve(rPoints);
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


Eigen::VectorXi NuTo::NURBSCurve::GetElementControlPointIDsGlobal(int rElementID, const Eigen::MatrixXi &rNodeIDs) const
{
    assert((rElementID >= 0) && (rElementID < GetNumIGAElements()));

    int knotID = GetElementFirstKnotID(rElementID);

    Eigen::VectorXi ids(mDegree+1);

    for(int i = 0; i < mDegree + 1; i++)
    {
        ids(i) = rNodeIDs(knotID - mDegree + i);
    }

    return ids;
}

Eigen::VectorXi NuTo::NURBSCurve::GetParameterControlPoints(double rParameter)
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    Eigen::VectorXi controlPointIDs(mDegree+1);
    for(int j = 0; j <= mDegree; j++) controlPointIDs(j) = spanIdx-mDegree+j;

    return controlPointIDs;
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

int NuTo::NURBSCurve::GetBasisFunctionIndexNotVanishing(double rParameter) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    return spanIdx-mDegree;
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
    rParameters.setZero(numPoints);

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

Eigen::VectorXd NuTo::NURBSCurve::CurvePoint(double rParameter, int rDerivative) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rParameter, mDegree, mKnots);
    Eigen::VectorXd basisFunctions = ShapeFunctionsIGA::BasisFunctionsAndDerivativesRat(rDerivative, rParameter, spanIdx, mDegree, mKnots, mWeights);

    Eigen::VectorXd coordinates(GetDimension());
    coordinates.setZero(GetDimension());

    for(int i = 0; i < GetDimension(); i++)
        for(int j = 0; j <= mDegree; j++)
            coordinates(i) += basisFunctions(j)*mControlPoints(spanIdx-mDegree+j, i);

    return coordinates;
}


Eigen::MatrixXd NuTo::NURBSCurve::CurvePoints(const Eigen::VectorXd &rParameter, int rDerivative) const
{
    int numParameters = rParameter.rows();
    Eigen::MatrixXd result(numParameters, GetDimension());
    for(int i = 0; i < numParameters; i++)
    {
        result.row(i) = CurvePoint(rParameter[i], rDerivative).transpose();
    }

    return result;
}

void NuTo::NURBSCurve::InsertKnot(double rKnotToInsert, int rMultiplicity)
{
    int k = ShapeFunctionsIGA::FindSpan(rKnotToInsert, mDegree, mKnots);
    int initialMultiplicity = GetMultiplicityOfKnot(rKnotToInsert);

    assert(initialMultiplicity + rMultiplicity <= 1);
    //assert(initialMultiplicity + rMultiplicity <= mDegree);

    // new knot vector
    Eigen::VectorXd newKnots(GetNumKnots() + rMultiplicity);

    for(int i = 0 ; i <= k; i++)               newKnots(i) = mKnots(i);
    for(int i = 1 ; i <= rMultiplicity; i++)   newKnots(k + i) = rKnotToInsert;
    for(int i = k + 1; i < GetNumKnots(); i++) newKnots(rMultiplicity + i) = mKnots(i);

    // since we're dealing with NURBS curves, a projection is needed ...
    Eigen::MatrixXd controlPointsProjected(mControlPoints.rows(), GetDimension()+1);
    controlPointsProjected.block(0, 0, mControlPoints.rows(), mControlPoints.cols()) = mControlPoints;

    int col = mControlPoints.cols();
    for(int i = 0; i < mControlPoints.rows(); i++)
        controlPointsProjected(i, col) = mWeights(i);

    for(int i = 0; i < controlPointsProjected.rows(); i++)
        for(int j = 0; j < controlPointsProjected.cols() - 1; j++)
            controlPointsProjected(i, j) *= mWeights(i);

    Eigen::MatrixXd newControlPoints(GetNumControlPoints() + rMultiplicity, GetDimension()+1);
    Eigen::MatrixXd Rw(mDegree+1, GetDimension() + 1);

    for(int i = 0; i <= k - mDegree; i++) newControlPoints.row(i) = controlPointsProjected.row(i);
    for(int i = k - initialMultiplicity; i < GetNumControlPoints(); i++) newControlPoints.row(i + rMultiplicity) = controlPointsProjected.row(i);
    for(int i = 0; i <= mDegree - initialMultiplicity; i++) Rw.row(i) = controlPointsProjected.row(k - mDegree + i);

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

    mWeights = newControlPoints.col(newControlPoints.cols() - 1);
    mControlPoints.resize(newControlPoints.rows(), GetDimension());
    for(int i = 0; i < newControlPoints.rows(); i++) mControlPoints.row(i) = newControlPoints.block(i,0, 1, GetDimension())/mWeights(i);
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

    for(int i = 0; i < controlPointsProjected.rows(); i++)
        for(int j = 0; j < controlPointsProjected.cols() - 1; j++)
            controlPointsProjected(i, j) *= mWeights(i);

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

void NuTo::NURBSCurve::findMinimalDistance(const Eigen::VectorXd &rCoordinatesSlave, double &rParameterStartMaster) const
{
    double tol = 1.e-10;
    double error = 1.;
    int maxNumIter = 100;
    int numIter = 0;

    if((rCoordinatesSlave - mControlPoints.block(GetNumControlPoints()-1, 0, 1, 2).transpose()).norm()< 1.e-10 )
        rParameterStartMaster = 1.;
    if((rCoordinatesSlave - mControlPoints.block(0, 0, 1, 2).transpose()).norm()< 1.e-10 )
        rParameterStartMaster = 0.;
    while(error > tol && numIter < maxNumIter)
    {
        // ==> function (dprime)
        Eigen::VectorXd coordinatesMaster = CurvePoint(rParameterStartMaster, 0);
        Eigen::VectorXd r = rCoordinatesSlave - coordinatesMaster;
        //if(r.norm() < tol) break;
        double b = 0.;
        Eigen::VectorXd fistDer   = CurvePoint(rParameterStartMaster, 1);
        Eigen::VectorXd secondDer = CurvePoint(rParameterStartMaster, 2);

        b = r.dot(fistDer);

        // ==> derivative
        double A = r.dot(secondDer) - fistDer.dot(fistDer);

        // ==> iteration step
        double incr = -b/A;
        rParameterStartMaster += incr;

        error = std::fabs(b);
        numIter++;
    }
//    std::cout << "Number of iterations: " << numIter << std::endl;
}

Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> NuTo::NURBSCurve::buildIGAStructure(NuTo::Structure &rStructure,
                                                                                                       const std::set<NuTo::Node::eDof> &rSetOfDOFS,
                                                                                                       int rGroupElements,
                                                                                                       int rGroupNodes,
                                                                                                       const std::string &rInterpolation,
                                                                                                       Eigen::VectorXi &nodeIDs) const
{
    nodeIDs.resize(GetNumControlPoints());

    // create dofs
    for(int i = 0; i < GetNumControlPoints(); i++)
    {
        int id = rStructure.NodeCreateDOFs(rSetOfDOFS, GetControlPoint(i));
        nodeIDs(i) = id;
        rStructure.GroupAddNode(rGroupNodes, id);
//        std::cout << "Node (id, coord): " <<  id << ", " << GetControlPoint(i).transpose() << std::endl;
    }

//    std::cout << "nodeIDs: " << nodeIDs.transpose() << std::endl;

    Eigen::VectorXi degree(1);
    degree << mDegree;

    int rNewInterpolation = rStructure.InterpolationTypeCreate(rInterpolation);
    std::vector<Eigen::VectorXd> vecKnots;
    vecKnots.push_back(GetKnotVector());
    for(auto &it : rSetOfDOFS) rStructure.InterpolationTypeAdd(rNewInterpolation, it, NuTo::Interpolation::eTypeOrder::SPLINE, degree, vecKnots, mWeights);

    int numElements = GetNumIGAElements();
    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elements(1, numElements);
    Eigen::VectorXi elementIncidence(mDegree + 1);
    for(int element = 0; element < numElements; element++)
    {
        elementIncidence = GetElementControlPointIDsGlobal(element, nodeIDs);
//        std::cout << "elementIncidence: " <<  elementIncidence.transpose() << std::endl;
        Eigen::MatrixXd knots = GetElementKnots(element);
        Eigen::VectorXi knotIds = GetElementKnotIDs(element);
        int id = rStructure.ElementCreate(rNewInterpolation, elementIncidence, knots, knotIds);
        rStructure.GroupAddElement(rGroupElements, id);
        elements(0, element) = std::pair<int, int>(id, -1);
    }

    return elements;
}


