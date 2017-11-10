#include "nuto/mechanics/IGA/NURBSSurface.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/math/FullMatrix.h"

NuTo::NURBSSurface::NURBSSurface(const Eigen::Vector2i& rDegree, const Eigen::VectorXd& rKnotsX,
                                 const Eigen::VectorXd& rKnotsY,
                                 const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>& rControlPoints,
                                 const Eigen::MatrixXd& rWeights)

    : mKnotsX(rKnotsX)
    , mKnotsY(rKnotsY)
    , mDegree(rDegree)
    , mWeights(rWeights)
    , mControlPoints(rControlPoints)
{
}

NuTo::NURBSSurface::NURBSSurface(const Eigen::Vector2i& rDegree,
                                 const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                                 Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& AInv)

    : mDegree(rDegree)
{
    // TODO
}

void NuTo::NURBSSurface::ParametrizationChordLengthMethod(const Eigen::MatrixXd& rPoints, Eigen::VectorXd& rParametersX,
                                                          Eigen::VectorXd& rParametersY)
{
    // TODO
}

int NuTo::NURBSSurface::GetNumControlPoints(int dir) const
{
    assert(dir == 0 || dir == 1);
    return ((dir == 0) ? mControlPoints.cols() : mControlPoints.rows());
}

int NuTo::NURBSSurface::GetNumControlPoints() const
{
    return mControlPoints.rows() * mControlPoints.cols();
}

Eigen::VectorXd NuTo::NURBSSurface::GetControlPoint(int rControlPointIDY, int rControlPointIDX) const
{
    assert(rControlPointIDY >= 0 && rControlPointIDY < GetNumControlPoints(1));
    assert(rControlPointIDX >= 0 && rControlPointIDX < GetNumControlPoints(0));
    return mControlPoints(rControlPointIDY, rControlPointIDX);
}

int NuTo::NURBSSurface::GetNumKnots(int dir) const
{
    assert(dir == 0 || dir == 1);
    return ((dir == 0) ? mKnotsX.rows() : mKnotsY.rows());
}
const Eigen::VectorXd& NuTo::NURBSSurface::GetKnotVector(int dir) const
{
    assert(dir == 0 || dir == 1);
    return ((dir == 0) ? mKnotsX : mKnotsY);
}

int NuTo::NURBSSurface::GetElementFirstKnotID(int rElementIDinDir, int dir) const
{
    assert(dir == 0 || dir == 1);

    int elementID = 0;
    int knotID = -1;

    const Eigen::VectorXd& knots = (dir == 0) ? mKnotsX : mKnotsY;

    for (int i = 1; i < knots.rows(); i++)
    {
        if (knots(i - 1) < knots(i))
            elementID++;
        if (rElementIDinDir == elementID - 1)
        {
            knotID = i - 1;
            break;
        }
    }
    if (knotID == -1)
        throw Exception("[NURBSSurface::GetElementFirstKnotID] ElementId not found.");

    return knotID;
}

Eigen::MatrixXd NuTo::NURBSSurface::GetElementKnots(int rElementIDX, int rElementIDY) const
{
    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);
    Eigen::MatrixXd knots(2, 2);
    knots << mKnotsX(knotIDX), mKnotsX(knotIDX + 1), mKnotsY(knotIDY), mKnotsY(knotIDY + 1);

    return knots;
}

Eigen::VectorXi NuTo::NURBSSurface::GetElementKnotIDs(int rElementIDX, int rElementIDY) const
{
    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);
    Eigen::VectorXi knots(2);
    knots << knotIDX, knotIDY;

    return knots;
}

int NuTo::NURBSSurface::GetMultiplicityOfKnot(double rKnot, int dir) const
{
    const Eigen::VectorXd& knots = (dir == 0) ? mKnotsX : mKnotsY;
    int spanIdx = ShapeFunctionsIGA::FindSpan(rKnot, mDegree(dir), knots);

    int count = 0;
    for (int i = spanIdx; i > 0; i--)
        if (knots(i) == rKnot)
            count++;
        else
            break;

    return count;
}

int NuTo::NURBSSurface::GetNumIGAElements(int dir) const
{
    assert(dir == 0 || dir == 1);
    const Eigen::VectorXd& knots = (dir == 0) ? mKnotsX : mKnotsY;

    int numElements = 0;

    for (int i = 1; i < knots.rows(); i++)
    {
        if (knots(i - 1) < knots(i))
            numElements++;
    }

    return numElements;
}

int NuTo::NURBSSurface::GetNumIGAElements() const
{
    return GetNumIGAElements(0) * GetNumIGAElements(1);
}


Eigen::VectorXi NuTo::NURBSSurface::GetElementControlPointIDs(int rElementIDX, int rElementIDY) const
{
    Eigen::VectorXi ids((mDegree(0) + 1) * (mDegree(1) + 1));

    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);

    int count = 0;
    for (int i = 0; i <= mDegree(1); i++)
    {
        int idy = knotIDY - mDegree(1) + i;
        for (int j = 0; j <= mDegree(0); j++)
        {
            int idx = knotIDX - mDegree(0) + j;
            ids(count) = idy * GetNumControlPoints(0) + idx;
            count++;
        }
    }

    return ids;
}

Eigen::VectorXi NuTo::NURBSSurface::GetElementControlPointIDsGlobal(int rElementIDX, int rElementIDY,
                                                                    const Eigen::MatrixXi& rNodeIDs) const
{
    Eigen::VectorXi ids((mDegree(0) + 1) * (mDegree(1) + 1));

    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);

    int count = 0;
    for (int i = 0; i <= mDegree(1); i++)
    {
        int idy = knotIDY - mDegree(1) + i;
        for (int j = 0; j <= mDegree(0); j++)
        {
            int idx = knotIDX - mDegree(0) + j;
            ids(count) = rNodeIDs(idy, idx);
            count++;
        }
    }

    return ids;
}

Eigen::MatrixXd NuTo::NURBSSurface::GetKnotIDControlPoints(const Eigen::Vector2i& rKnotIDs) const
{
    Eigen::MatrixXd coords((mDegree(0) + 1) * (mDegree(1) + 1), GetDimension());

    int count = 0;
    for (int i = 0; i <= mDegree(1); i++)
    {
        int idy = rKnotIDs(1) - mDegree(1) + i;
        for (int j = 0; j <= mDegree(0); j++)
        {
            int idx = rKnotIDs(0) - mDegree(0) + j;
            coords.row(count) = mControlPoints(idy, idx);
            count++;
        }
    }

    return coords;
}

Eigen::VectorXd NuTo::NURBSSurface::SurfacePoint(const Eigen::Vector2d& rParameter) const
{
    Eigen::Vector2i span;
    span(0) = ShapeFunctionsIGA::FindSpan(rParameter(0), mDegree(0), mKnotsX);
    span(1) = ShapeFunctionsIGA::FindSpan(rParameter(1), mDegree(1), mKnotsY);

    Eigen::VectorXd basisFunctions =
            ShapeFunctionsIGA::BasisFunctions2DRat(rParameter, span, mDegree, mKnotsX, mKnotsY, mWeights);

    Eigen::VectorXd coordinates(GetDimension());

    Eigen::MatrixXd cpCoords = GetKnotIDControlPoints(span);

    for (int i = 0; i < GetDimension(); i++)
        coordinates(i) = basisFunctions.dot(cpCoords.col(i));

    return coordinates;
}

Eigen::MatrixXd NuTo::NURBSSurface::SurfacePoints(const Eigen::MatrixXd& rParameter) const
{
    assert(rParameter.cols() == 2);

    Eigen::MatrixXd coordinates(rParameter.rows(), GetDimension());

    for (int i = 0; i < rParameter.rows(); i++)
    {
        coordinates.row(i) = SurfacePoint(rParameter.row(i));
    }

    return coordinates;
}

void NuTo::NURBSSurface::InsertKnot(double rKnotToInsert, int rMultiplicity, int dir)
{
    assert(dir == 0 || dir == 1);

    Eigen::VectorXd& knots = (dir == 0) ? mKnotsX : mKnotsY;

    int k = ShapeFunctionsIGA::FindSpan(rKnotToInsert, mDegree(dir), knots);
    int initialMultiplicity = GetMultiplicityOfKnot(rKnotToInsert, dir);

    assert(initialMultiplicity + rMultiplicity <= mDegree(dir));


    // new knot vector
    Eigen::VectorXd newKnots(GetNumKnots(dir) + rMultiplicity);

    for (int i = 0; i <= k; i++)
        newKnots(i) = knots(i);
    for (int i = 1; i <= rMultiplicity; i++)
        newKnots(k + i) = rKnotToInsert;
    for (int i = k + 1; i < GetNumKnots(dir); i++)
        newKnots(rMultiplicity + i) = knots(i);

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> newControlPoints;
    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> Rw(mDegree(dir) + 1, 1);
    Eigen::MatrixXd alpha(mDegree(dir) - initialMultiplicity, rMultiplicity + 1);

    if (dir == 0)
    {
        // new control points
        newControlPoints.resize(GetNumControlPoints(1), GetNumControlPoints(dir) + rMultiplicity);

        int L = 0;
        for (int j = 1; j <= rMultiplicity; j++)
        {
            L = k - mDegree(dir) + j;
            for (int i = 0; i <= mDegree(dir) - j - initialMultiplicity; i++)
                alpha(i, j) = (rKnotToInsert - knots(L + i)) / (knots(i + k + 1) - knots(L + i));
        }

        for (int row = 0; row < GetNumControlPoints(1); row++)
        {
            for (int i = 0; i <= k - mDegree(dir); i++)
                newControlPoints(row, i) = mControlPoints(row, i);
            for (int i = k - initialMultiplicity; i < GetNumControlPoints(dir); i++)
                newControlPoints(row, i + rMultiplicity) = mControlPoints(row, i);
            for (int i = 0; i <= mDegree(dir) - initialMultiplicity; i++)
                Rw(i) = mControlPoints(row, k - mDegree(dir) + i);

            for (int j = 1; j <= rMultiplicity; j++)
            {
                L = k - mDegree(dir) + j;
                for (int i = 0; i <= mDegree(dir) - j - initialMultiplicity; i++)
                    Rw(i) = alpha(i, j) * Rw(i + 1) + (1.0 - alpha(i, j)) * Rw(i);
                newControlPoints(row, L) = Rw(0);
                newControlPoints(row, k + rMultiplicity - j - initialMultiplicity) =
                        Rw(mDegree(dir) - j - initialMultiplicity);
            }

            for (int i = L + 1; i < k - initialMultiplicity; i++)
                newControlPoints(row, i) = Rw(i - L);
        }
    }

    knots = newKnots;
    mControlPoints = newControlPoints;
}

void NuTo::NURBSSurface::RefineKnots(const Eigen::VectorXd& rKnotsToInsert, int dir)
{
    assert(dir == 0 || dir == 1);

    Eigen::VectorXd& knots = (dir == 0) ? mKnotsX : mKnotsY;

    int numInsert = rKnotsToInsert.rows();
    int begin = ShapeFunctionsIGA::FindSpan(rKnotsToInsert(0), mDegree(dir), knots);
    int end = ShapeFunctionsIGA::FindSpan(rKnotsToInsert(numInsert - 1), mDegree(dir), knots);
    end++;

    int dimension = GetDimension();
    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPointsProjected(mControlPoints.rows(),
                                                                                          mControlPoints.cols());

    Eigen::VectorXd temp(dimension + 1);
    for (int i = 0; i < mControlPoints.rows(); i++)
    {
        for (int j = 0; j < mControlPoints.cols(); j++)
        {
            temp.block(0, 0, dimension, 1) = mControlPoints(i, j);
            temp *= mWeights(i, j);
            temp(dimension) = mWeights(i, j);
            controlPointsProjected(i, j) = temp;
        }
    }

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> newControlPoints;

    // new knot vector
    Eigen::VectorXd newKnots(GetNumKnots(dir) + numInsert);

    for (int i = 0; i <= begin; i++)
        newKnots(i) = knots(i);
    for (int i = end + mDegree(dir); i < GetNumKnots(dir); i++)
        newKnots(numInsert + i) = knots(i);

    if (dir == 0)
    {
        // new control points
        newControlPoints.resize(GetNumControlPoints(1), GetNumControlPoints(dir) + numInsert);

        for (int row = 0; row < GetNumControlPoints(1); row++)
        {
            for (int j = 0; j <= begin - mDegree(dir); j++)
                newControlPoints(row, j) = controlPointsProjected(row, j);
            for (int j = end - 1; j < GetNumControlPoints(dir); j++)
                newControlPoints(row, j + numInsert) = controlPointsProjected(row, j);
        }

        int i = end + mDegree(dir) - 1;
        int k = end + mDegree(dir) + numInsert - 1;

        for (int j = numInsert - 1; j >= 0; j--)
        {
            while (rKnotsToInsert(j) <= knots(i) && i > begin)
            {
                newKnots(k) = knots(i);
                for (int row = 0; row < GetNumControlPoints(1); row++)
                    newControlPoints(row, k - mDegree(dir) - 1) = controlPointsProjected(row, i - mDegree(dir) - 1);
                k--;
                i--;
            }
            for (int row = 0; row < GetNumControlPoints(1); row++)
                newControlPoints(row, k - mDegree(dir) - 1) = newControlPoints(row, k - mDegree(dir));
            for (int l = 1; l <= mDegree(dir); l++)
            {
                int ind = k - mDegree(dir) + l;
                double alpha = newKnots(k + l) - rKnotsToInsert(j);
                if (std::fabs(alpha) == 0.0)
                {
                    for (int row = 0; row < GetNumControlPoints(1); row++)
                        newControlPoints(row, ind - 1) = newControlPoints(row, ind);
                }
                else
                {
                    alpha /= (newKnots(k + l) - knots(i - mDegree(dir) + l));
                    for (int row = 0; row < GetNumControlPoints(1); row++)
                        newControlPoints(row, ind - 1) =
                                alpha * newControlPoints(row, ind - 1) + (1.0 - alpha) * newControlPoints(row, ind);
                }
            }
            newKnots(k) = rKnotsToInsert(j);
            k--;
        }
    }
    if (dir == 1)
    {
        // new control points
        newControlPoints.resize(GetNumControlPoints(1) + numInsert, GetNumControlPoints(0));

        for (int col = 0; col < GetNumControlPoints(0); col++)
        {
            for (int j = 0; j <= begin - mDegree(dir); j++)
                newControlPoints(j, col) = controlPointsProjected(j, col);
            for (int j = end - 1; j < GetNumControlPoints(dir); j++)
                newControlPoints(j + numInsert, col) = controlPointsProjected(j, col);
        }

        int i = end + mDegree(dir) - 1;
        int k = end + mDegree(dir) + numInsert - 1;

        for (int j = numInsert - 1; j >= 0; j--)
        {
            while (rKnotsToInsert(j) <= knots(i) && i > begin)
            {
                newKnots(k) = knots(i);
                for (int col = 0; col < GetNumControlPoints(0); col++)
                    newControlPoints(k - mDegree(dir) - 1, col) = controlPointsProjected(i - mDegree(dir) - 1, col);
                k = k - 1;
                i = i - 1;
            }
            for (int col = 0; col < GetNumControlPoints(0); col++)
                newControlPoints(k - mDegree(dir) - 1, col) = newControlPoints(k - mDegree(dir), col);
            for (int l = 1; l <= mDegree(dir); l++)
            {
                int ind = k - mDegree(dir) + l;
                double alpha = newKnots(k + l) - rKnotsToInsert(j);
                if (std::fabs(alpha) == 0.0)
                {
                    for (int col = 0; col < GetNumControlPoints(0); col++)
                        newControlPoints(ind - 1, col) = newControlPoints(ind, col);
                }
                else
                {
                    alpha /= (newKnots(k + l) - knots(i - mDegree(dir) + l));
                    for (int col = 0; col < GetNumControlPoints(0); col++)
                        newControlPoints(ind - 1, col) =
                                alpha * newControlPoints(ind - 1, col) + (1.0 - alpha) * newControlPoints(ind, col);
                }
            }
            newKnots(k) = rKnotsToInsert(j);
            k--;
        }
    }
    knots = newKnots;

    mControlPoints.resize(newControlPoints.rows(), newControlPoints.cols());
    mWeights.resize(newControlPoints.rows(), newControlPoints.cols());
    temp(dimension);
    for (int i = 0; i < newControlPoints.rows(); i++)
    {
        for (int j = 0; j < newControlPoints.cols(); j++)
        {
            temp = newControlPoints(i, j).block(0, 0, dimension, 1);
            temp /= newControlPoints(i, j)(dimension);
            mControlPoints(i, j) = temp;
            mWeights(i, j) = newControlPoints(i, j)(dimension);
        }
    }
}

void NuTo::NURBSSurface::DuplicateKnots(int dir)
{
    assert(dir == 0 || dir == 1);

    Eigen::VectorXd& knots = (dir == 0) ? mKnotsX : mKnotsY;

    Eigen::VectorXd knotsToInsert(GetNumIGAElements(dir));
    for (int i = 0; i < knotsToInsert.rows(); i++)
    {
        int beginID = GetElementFirstKnotID(i, dir);
        knotsToInsert(i) = (knots(beginID + 1) + knots(beginID)) / 2.;
    }
    RefineKnots(knotsToInsert, dir);
}

Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic>
NuTo::NURBSSurface::buildIGAStructure(NuTo::Structure& rStructure, const std::set<NuTo::Node::eDof>& rSetOfDOFS,
                                      int rGroupElements, int rGroupNodes, const std::string& rInterpolation)
{
    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elements(GetNumIGAElements(1),
                                                                                GetNumIGAElements(0));

    Eigen::MatrixXi nodeIDs(GetNumControlPoints(1), GetNumControlPoints(0));

    // create dofs
    for (int i = 0; i < GetNumControlPoints(1); i++)
        for (int j = 0; j < GetNumControlPoints(0); j++)
        {
            int id = rStructure.NodeCreateDOFs(rSetOfDOFS, GetControlPoint(i, j));
            nodeIDs(i, j) = id;
            rStructure.GroupAddNode(rGroupNodes, id);
        }

    int rNewInterpolation = rStructure.InterpolationTypeCreate(rInterpolation);
    std::vector<Eigen::VectorXd> vecKnots;
    vecKnots.push_back(GetKnotVector(0));
    vecKnots.push_back(GetKnotVector(1));
    for (auto& it : rSetOfDOFS)
        rStructure.InterpolationTypeAdd(rNewInterpolation, it, NuTo::Interpolation::eTypeOrder::SPLINE, mDegree,
                                        vecKnots, mWeights);


    Eigen::VectorXi elementIncidence((mDegree(0) + 1) * (mDegree(1) + 1));
    for (int elementY = 0; elementY < GetNumIGAElements(1); elementY++)
    {
        for (int elementX = 0; elementX < GetNumIGAElements(0); elementX++)
        {
            elementIncidence = GetElementControlPointIDsGlobal(elementX, elementY, nodeIDs);
            int id = rStructure.ElementCreate(rNewInterpolation, elementIncidence, GetElementKnots(elementX, elementY),
                                              GetElementKnotIDs(elementX, elementY));
            rStructure.GroupAddElement(rGroupElements, id);
            elements(elementY, elementX) = std::pair<int, int>(id, -1);
        }
    }

    return elements;
}
