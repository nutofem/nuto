#include "BoostUnitTest.h"

#include "mechanics/timeIntegration/TimeControl.h"

#include <iostream>
#include "base/Exception.h"

using namespace NuTo;

// %%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool CheckExceptionMessageCorrect(Exception const& ex, std::string expectedMsg)
{
    std::string refinedErrorMessage = ex.ErrorMessage().substr(ex.ErrorMessage().find_first_of("]") + 2);

    if (expectedMsg.compare(refinedErrorMessage))
    {
        std::cout << std::endl
                  << "WRONG EXCEPTION TRIGGERED!" << std::endl
                  << "Expected exception message:" << std::endl
                  << expectedMsg << std::endl
                  << "Triggered exception message:" << std::endl
                  << refinedErrorMessage << std::endl
                  << std::endl;
        return false;
    }
    return true;
}


// %%% Tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BOOST_AUTO_TEST_CASE(Debug_Fuctionality)
{
    TimeControl timeControl;

    std::string expectedMsg = "";
    std::function<bool(Exception const&)> CheckException = [&expectedMsg](Exception const& ex) -> bool {
        return CheckExceptionMessageCorrect(ex, expectedMsg);
    };

    expectedMsg = "Current timestep is 0 or negative!";
    timeControl.SetTimeStepFunction([](TimeControl&, double, double, bool) -> double { return 0.0; });
    BOOST_CHECK_EXCEPTION(timeControl.Proceed(), Exception, CheckException);
    timeControl.SetTimeStepFunction([](TimeControl&, double, double, bool) -> double { return -1.0; });
    BOOST_CHECK_EXCEPTION(timeControl.Proceed(), Exception, CheckException);

    expectedMsg = "Timestep must be larger than 0!";
    BOOST_CHECK_EXCEPTION(timeControl.SetTimeStep(0.0), Exception, CheckException);
    BOOST_CHECK_EXCEPTION(timeControl.SetTimeStep(-1.0), Exception, CheckException);

    expectedMsg = "Maximal timestep must be a positive number!";
    BOOST_CHECK_EXCEPTION(timeControl.SetMaxTimeStep(0.0), Exception, CheckException);
    BOOST_CHECK_EXCEPTION(timeControl.SetMaxTimeStep(-1.0), Exception, CheckException);

    timeControl.SetMinTimeStep(1.0);
    timeControl.SetMaxTimeStep(5.0);
    expectedMsg = "Minimal timestep must be smaller than maximal Timestep!";
    BOOST_CHECK_EXCEPTION(timeControl.SetMinTimeStep(6.0), Exception, CheckException);
    expectedMsg = "Maximal timestep must be bigger than minimal Timestep!";
    BOOST_CHECK_EXCEPTION(timeControl.SetMaxTimeStep(0.5), Exception, CheckException);


    timeControl.SetTimeStepFunction(
            [](TimeControl& timeControl, double, double, bool) -> double { return timeControl.GetTimeStep(); });
    timeControl.SetTimeStep(1.0);
    timeControl.Proceed();
    timeControl.Proceed();

    expectedMsg = "Final time must be larger than current time!";
    BOOST_CHECK_EXCEPTION(timeControl.SetTimeFinal(1.0), Exception, CheckException);


    expectedMsg = "No convergence with the current maximum number of "
                  "iterations, either use automatic time stepping, "
                  "reduce the time step or the minimal line search cut "
                  "back factor. In case you provided a custom timestepping function "
                  "and intend to use some kind of automatic timestepping, call "
                  "the RestorePreviosTime() function of the time control before reducing "
                  "the timestep.";
    BOOST_CHECK_EXCEPTION(timeControl.AdjustTimestep(10, 10, false), Exception, CheckException);
}


BOOST_AUTO_TEST_CASE(Timestepping)
{
    // Equidistant timestepping
    {
        TimeControl timeControl;
        timeControl.SetTimeFinal(20);


        timeControl.SetTimeStep(10.);
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 10.0);
        BOOST_CHECK_EQUAL(timeControl.Finished(), false);

        timeControl.SetTimeStep(5.0);
        timeControl.UseEquidistantTimestepping();

        timeControl.Proceed();
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 20.0);
        BOOST_CHECK_EQUAL(timeControl.Finished(), true);
    }


    // Automatic timestepping
    {
        TimeControl timeControl;
        timeControl.SetTimeFinal(6);

        timeControl.SetMinTimeStep(1.0);
        timeControl.SetMaxTimeStep(4.0);
        timeControl.SetTimeStep(2.0);

        timeControl.UseDefaultAutomaticTimestepping();

        // timestep should be adjusted to 3.0
        timeControl.AdjustTimestep(1, 20, true);
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 3.0);
        BOOST_CHECK_EQUAL(timeControl.Finished(), false);
        BOOST_CHECK_EQUAL(timeControl.GetTimeStep(), 3.0);

        timeControl.SetTimeStep(2.0);
        // timestep should remain 2.0
        timeControl.AdjustTimestep(10, 20, true);
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetTimeStep(), 2.0);
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 5.0);
        BOOST_CHECK_EQUAL(timeControl.Finished(), false);

        // Proceed to time = 7.0
        timeControl.Proceed();
        // No convergence -> timestep should be adjusted to 1.0 and current time should be set back to 5.0
        timeControl.AdjustTimestep(1, 20, false);
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetTimeStep(), 1.0);
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 6.0);
        BOOST_CHECK_EQUAL(timeControl.Finished(), true);
    }

    // Automatic timestepping, min time step
    {
        TimeControl timeControl;
        timeControl.SetMinTimeStep(1.0);
        timeControl.SetTimeFinal(6);
        timeControl.SetTimeStep(1.5);

        timeControl.UseDefaultAutomaticTimestepping();

        // should decrease the time step below min time step
        BOOST_CHECK_THROW(timeControl.AdjustTimestep(1, 20, false), NuTo::Exception);
    }

    // floating point accuracy of TimeControl::Finished()
    {
        TimeControl timeControl;
        constexpr int n = 100000;
        timeControl.SetTimeFinal(n);
        timeControl.SetTimeStep(1);
        for (int i = 0; i < n; ++i)
            timeControl.Proceed();
        BOOST_CHECK(timeControl.Finished());
    }
}
