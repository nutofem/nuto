#ifndef TIMEDEPENDENCYMATRIX
#define TIMEDEPENDENCYMATRIX


#include <nuto/math/FullMatrix.h>

namespace NuTo
{
//! @author Volker Hirthammer
//! @date July 09, 2015
//! @brief ...
class TimeDependencyMatrix : public TimeDependencyBase
{
public:
    TimeDependencyMatrix(const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rTimeDependencyMatrix)
        : TimeDependencyBase(),
          mTimeDependencyMatrix(rTimeDependencyMatrix)
    {}

    virtual ~TimeDependencyMatrix(){}

    virtual double GetTimeDependentFactor(double rTime) override
    {
        //calculate the two corresponding time steps between which a linear interpolation is performed
        if (mTimeDependencyMatrix.GetNumRows()!=0)
        {
            int curStep(0);
            while (mTimeDependencyMatrix(curStep,0)<rTime && curStep<mTimeDependencyMatrix.GetNumRows()-1)
                curStep++;
            if (curStep==0)
                curStep++;

            //extract the two data points
            double s1 = mTimeDependencyMatrix(curStep-1,1);
            double s2 = mTimeDependencyMatrix(curStep,1);
            double t1 = mTimeDependencyMatrix(curStep-1,0);
            double t2 = mTimeDependencyMatrix(curStep,0);

            return s1 + (s2-s1)/(t2-t1) * (rTime-t1);
        }
        return 0;

    }

private:
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mTimeDependencyMatrix;
};
} // namespace NuTo

#endif // TIMEDEPENDENCYMATRIX

