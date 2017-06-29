#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#include <boost/ptr_container/serialize_ptr_list.hpp>
#endif // ENABLE_SERIALIZATION

#include "optimize/NonlinearSolverBase.h"
#include "base/Exception.h"

using namespace std;

NuTo::NonlinearSolverBase::NonlinearSolverBase()
{
    mTolResidual = 1.0e-12;
    mTolSolution = numeric_limits<double>::epsilon();
    mMaxIterationsNumber = 100;
    mResidualFunction = nullptr;
    mAssignResidual = false;
}

Eigen::MatrixXd NuTo::NonlinearSolverBase::DResidualNum(Eigen::VectorXd rUnknown, Eigen::VectorXd& rFvec) const
{
    const double EPS = 1.0e-8;
    int n = rUnknown.rows();
    Eigen::MatrixXd deriv(n, n);
    Eigen::VectorXd xh = rUnknown;

    if (mResidualFunction == nullptr && mAssignResidual == false)
    {
        throw Exception(__PRETTY_FUNCTION__, "The pointer to the residual function is required.");
    }

    for (int j = 0; j < n; j++)
    {
        double temp = xh[j];
        double h = EPS * std::abs(temp);
        if (h == 0.0)
            h = EPS;
        xh[j] = temp + h;
        h = xh[j] - temp;
        Eigen::VectorXd f = (mResidualFunctionBoost)(this->mParameter, xh);
        xh[j] = temp;
        for (int i = 0; i < n; i++)
            deriv(i, j) = (f[i] - rFvec[i]) / h;
    }
    return deriv;
}

double NuTo::NonlinearSolverBase::Fmin(Eigen::VectorXd rUnknown, Eigen::VectorXd& rFvec) const
{
    if (mResidualFunction == nullptr && mAssignResidual == false)
    {
        throw Exception(__PRETTY_FUNCTION__, "The pointer to the residual function is required.");
    }

    rFvec = (mResidualFunctionBoost)(this->mParameter, rUnknown);
    double sum = 0;
    sum = rFvec.dot(rFvec);
    return 0.5 * sum;
}

void NuTo::NonlinearSolverBase::Info() const {}

#ifdef ENABLE_SERIALIZATION
template void NuTo::NonlinearSolverBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NonlinearSolverBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of NonlinearSolverBase" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mTolResidual)
       & BOOST_SERIALIZATION_NVP(mTolSolution)
       & BOOST_SERIALIZATION_NVP(mMaxIterationsNumber)
       & BOOST_SERIALIZATION_NVP(mParameter);
//       & BOOST_SERIALIZATION_NVP(mResidualFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of NonlinearSolverBase" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION
