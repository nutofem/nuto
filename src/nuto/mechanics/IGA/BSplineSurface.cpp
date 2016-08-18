#include "nuto/mechanics/IGA/BSplineSurface.h"

NuTo::BSplineSurface::BSplineSurface(const Eigen::Vector2i &rDegree,
                                     const Eigen::VectorXd &rKnotsX,
                                     const Eigen::VectorXd &rKnotsY,
                                     const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> &rControlPoints)

    : mKnotsX(rKnotsX), mKnotsY(rKnotsY), mDegree(rDegree), mControlPoints(rControlPoints)
{}

NuTo::BSplineSurface::BSplineSurface(const Eigen::Vector2i &rDegree,
                                     const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> &rPoints,
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &AInv)

    : mDegree(rDegree)
{
// TODO
}

void NuTo::BSplineSurface::ParametrizationChordLengthMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                                                            FullVector<double, Eigen::Dynamic>& rParametersX,
                                                            FullVector<double, Eigen::Dynamic>& rParametersY)
{
// TODO
}

int NuTo::BSplineSurface::GetNumControlPoints(int dir) const
{
    assert(dir == 0 || dir == 1);
    return ((dir == 0) ? mControlPoints.rows() : mControlPoints.cols());
}

int NuTo::BSplineSurface::GetNumControlPoints() const
{
    return mControlPoints.rows()*mControlPoints.cols();
}

Eigen::VectorXd NuTo::BSplineSurface::GetControlPoint(int rControlPointIDX, int rControlPointIDY) const
{
    assert(rControlPointIDX >= 0 && rControlPointIDX < GetNumControlPoints(0));
    assert(rControlPointIDY >= 0 && rControlPointIDY < GetNumControlPoints(1));
    return mControlPoints(rControlPointIDX, rControlPointIDY);
}

const Eigen::VectorXd& NuTo::BSplineSurface::GetKnotVector(int dir) const
{
    assert(dir == 0 || dir == 1);
    return ((dir == 0) ? mKnotsX : mKnotsY);
}

int NuTo::BSplineSurface::GetElementFirstKnotID(int rElementIDinDir, int dir) const
{
    assert(dir == 0 || dir == 1);

    int elementID = 0;
    int knotID = -1;

    const Eigen::VectorXd &knots = (dir == 0) ? mKnotsX : mKnotsY;

    for(int i = 1; i < knots.rows(); i++)
    {
        if (knots(i-1) < knots(i)) elementID++;
        if (rElementIDinDir == elementID-1)
        {
            knotID = i-1;
            break;
        }
    }
    if (knotID == -1)
        throw Exception("[BSplineSurface::GetElementFirstKnotID] ElementId not found.");

    return knotID;
}

Eigen::MatrixXd NuTo::BSplineSurface::GetElementKnots(int rElementIDX, int rElementIDY) const
{
    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);
    Eigen::MatrixXd knots(2,2);
    knots << mKnotsX(knotIDX) , mKnotsX(knotIDX+1),
             mKnotsY(knotIDY) , mKnotsY(knotIDY+1);

    return knots;
}

Eigen::VectorXi NuTo::BSplineSurface::GetElementKnotIDs(int rElementIDX, int rElementIDY) const
{
    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);
    Eigen::VectorXi knots(2);
    knots << knotIDX, knotIDY;

    return knots;
}

int NuTo::BSplineSurface::GetNumIGAElements(int dir) const
{
    assert(dir == 0 || dir == 1);
    const Eigen::VectorXd &knots = (dir == 0) ? mKnotsX : mKnotsY;

    int numElements = 0;

    for(int i = 1; i < knots.rows(); i++)
    {
        if(knots(i-1) < knots(i)) numElements++;
    }

    return numElements;
}


Eigen::VectorXi NuTo::BSplineSurface::GetElementControlPointIDs(int rElementIDX, int rElementIDY) const
{
    Eigen::VectorXi ids((mDegree(0)+1)*(mDegree(1)+1));

    int knotIDX = GetElementFirstKnotID(rElementIDX, 0);
    int knotIDY = GetElementFirstKnotID(rElementIDY, 1);

    int count = 0;
    for(int i = 0; i <= mDegree(1); i++)
    {
        int idy = knotIDY - mDegree(1) + i;
        for(int j = 0; j <= mDegree(0); j++)
        {
            int idx = knotIDX - mDegree(0) + j;
            ids(count) = idy*GetNumControlPoints(0) + idx;
            count++;
        }
    }

    return ids;
}
