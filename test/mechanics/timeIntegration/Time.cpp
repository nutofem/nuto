#include "BoostUnitTest.h"

#include "mechanics/timeIntegration/Time.h"

#include <iostream>

using namespace NuTo;

// %%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool CheckExceptionMessageCorrect(MechanicsException const& ex,std::string expectedMsg)
{
    std::string refinedErrorMessage = ex.ErrorMessage().substr(ex.ErrorMessage().find_first_of("]")+2);

    if (expectedMsg.compare(refinedErrorMessage))
    {
        std::cout << std::endl
                  << "WRONG EXCEPTION TRIGGERED!" << std::endl
                  << "Expected exception message:" << std::endl << expectedMsg << std::endl
                  << "Triggered exception message:" << std::endl << refinedErrorMessage << std::endl << std::endl;
        return false;
    }
    return true;
}


// %%% Tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef NDEBUG
BOOST_AUTO_TEST_CASE(Debug_Only_Fuctionality)
{
    Time time;

    std::string expectedMsg = "";
    std::function<bool(MechanicsException const&)> CheckException
            = [&expectedMsg](MechanicsException const& ex)->bool
                {
                    return CheckExceptionMessageCorrect(ex,expectedMsg);
                };

    time.SetTimestepFunction([](double curTime)->double{return curTime;});

    expectedMsg = "The current proceed function of the class does not increase the current time!";
    BOOST_CHECK_EXCEPTION(time.Proceed(),MechanicsException,CheckException);



}
#endif





BOOST_AUTO_TEST_CASE(Equidistant_Timestepping)
{
    std::string expectedMsg = "";
    std::function<bool(MechanicsException const&)> CheckException
            = [&expectedMsg](MechanicsException const& ex)->bool
                {
                    return CheckExceptionMessageCorrect(ex,expectedMsg);
                };
    {
        Time time;

        expectedMsg = "Timestep must be a positive number!";
        BOOST_CHECK_EXCEPTION(time.SetEquidistantTimestepping(0),MechanicsException,CheckException);
        BOOST_CHECK_EXCEPTION(time.SetEquidistantTimestepping(-1337.0),MechanicsException,CheckException);

        time.SetEquidistantTimestepping(1);
        time.Proceed();
        time.Proceed();
        BOOST_CHECK_EQUAL(time.Proceed(), 3.0);

        time.SetEquidistantTimestepping(5.5);

        time.Proceed();
        time.Proceed();
        BOOST_CHECK_EQUAL(time.Proceed(), 19.5);
    }
    {
        Time time;

    }

}
