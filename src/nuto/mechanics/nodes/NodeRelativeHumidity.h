#ifndef NODERELATIVEHUMIDITY_H
#define NODERELATIVEHUMIDITY_H


#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{

//! @author Volker Hirthammer, BAM
//! @brief class for relative humidity dofs
//! @date May 2015
template<int TNumRelativeHumidity, int TNumTimeDerivatives>
class NodeRelativeHumidity: public virtual NodeBase
{
public:
    //! @brief constructor
    NodeRelativeHumidity()
    {
        static_assert(TNumRelativeHumidity >= 0 and TNumRelativeHumidity <= 1, 	"The node must have 0 or 1 relative humidity.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, 	"The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mRelativeHumidity[iTimeDerivative] = Eigen::Matrix<double, TNumRelativeHumidity, 1>::Zero();

        mDofRelativeHumidity.fill(0);
    }

    //! @brief copy constructor
    NodeRelativeHumidity(const NodeRelativeHumidity&) = default;

    //! @brief operator =
    NodeRelativeHumidity& operator=(const NodeRelativeHumidity&) = default;

    //! @brief destructor
    virtual ~NodeRelativeHumidity()
    {}

    //! @brief gets the number of relative humidity dofs
    int GetNumRelativeHumidity() const override
    {
        return TNumRelativeHumidity;
    }

    //! @brief gets the relative humidity DOF number
    int GetDofRelativeHumidity() const override
    {
        return mDofRelativeHumidity[0];
    }

    //! @brief gets the relative humidity value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    double GetRelativeHumidity(int rTimeDerivative) const override
    {
        assert(rTimeDerivative <= 2);
        return mRelativeHumidity[rTimeDerivative](0);
    }

    //! @brief sets the relative humidity value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    //! @param rRelativeHumidity ... relative humidity
    void SetRelativeHumidity(int rTimeDerivative, double rRelativeHumidity) override
    {
        assert(rTimeDerivative <= 2);
        mRelativeHumidity[rTimeDerivative](0) = rRelativeHumidity;
    }

protected:
    std::array<Eigen::Matrix<double, TNumRelativeHumidity, 1>, TNumTimeDerivatives + 1> mRelativeHumidity;
    std::array<int, TNumRelativeHumidity> mDofRelativeHumidity;




};

} // namespace NuTo;




#endif // NODERELATIVEHUMIDITY_H
