#pragma once
#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <memory>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <cmath>
#include <chrono>

namespace BenchmarkInternal
{

class Runner
{
public:
    Runner()
    {
        mWallTimes.reserve(1e7);
    }
    bool KeepRunningTime(double rSeconds)
    {
        mWallTimes.push_back(std::chrono::high_resolution_clock::now());
        return GetDurationInSeconds() < rSeconds;
    }

    bool KeepRunningIterations(size_t rNumIterations)
    {
        mWallTimes.push_back(std::chrono::high_resolution_clock::now());
        return mWallTimes.size() < rNumIterations + 1; // plus one since the first iteration is timed as well
    }

    std::vector<double> GetSeconds() const
    {
        if (mWallTimes.empty())
            return {0}; // case: runner not used

        std::vector<double> seconds(mWallTimes.size() - 1);
        for (size_t i = 0; i < seconds.size(); ++i)
            seconds[i] = std::chrono::duration<double>(mWallTimes[i + 1] - mWallTimes[i]).count();
        return seconds;
    }

private:
    double GetDurationInSeconds() const
    {
        return std::chrono::duration<double>(*mWallTimes.rbegin() - *mWallTimes.begin()).count();
    }
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> mWallTimes;
};


struct Result
{
    Result(std::string rTestCase, std::string rTestName, std::vector<double> rSeconds)
        : mTestCase(rTestCase)
        , mTestName(rTestName)
    {
        auto b = rSeconds.begin();
        auto e = rSeconds.end();

        mRuns = std::distance(b, e);
        mTimeTotal = std::accumulate(b, e, 0.);
        mTimeMean = mTimeTotal / mRuns;
        double variance = 0.;
        for (double second : rSeconds)
            variance += (second - mTimeMean) * (second - mTimeMean); // lambda...
        mTimeStdDev = std::sqrt(variance / mRuns);
    }

    std::string mTestCase;
    std::string mTestName;
    int mRuns;
    double mTimeTotal;
    double mTimeMean;
    double mTimeStdDev;
};


class Reporter
{
public:
    void AddTestResult(Result rTestResult)
    {
        std::string name = rTestResult.mTestCase + "::" + rTestResult.mTestName;
        mTestResults.push_back(rTestResult);
        std::cout << "Test " << std::setw(40) << std::left << name << " took " << rTestResult.mTimeTotal << "s"
                  << std::endl;
    }
    void Report()
    {
        WriteHeader();
        for (const auto& testResult : mTestResults)
            WriteTestResult(testResult);
        WriteFooter();
    }

protected:
    virtual void WriteHeader(){};
    virtual void WriteTestResult(const Result&){};
    virtual void WriteFooter(){};

private:
    std::vector<Result> mTestResults;
};

class ReporterCout : public Reporter
{
protected:
    virtual void WriteHeader() override
    {
        std::cout << "\n\n\n";
        std::cout << "|                                          |        [ns]       |     runs    | std dev[ns] |\n";
        std::cout << "|------------------------------------------|-------------------|-------------|-------------|\n";
    }
    virtual void WriteTestResult(const Result& rT) override
    {
        long NS = static_cast<long>(std::round(rT.mTimeMean * 1e9));
        long stdDevNS = static_cast<long>(std::round(rT.mTimeStdDev * 1e9));
        std::string nsFormatted = InsertThousand(NS, " ");
        std::string stdDevFormatted = InsertThousand(stdDevNS, " ");
        std::string testName = rT.mTestCase + "::" + rT.mTestName;
        std::cout << "| ";
        std::cout << std::left << std::setw(40) << testName << " |";
        std::cout << std::right << std::setw(18) << nsFormatted << " |";
        std::cout << std::right << std::setw(12) << rT.mRuns << " |";
        std::cout << std::right << std::setw(12) << stdDevFormatted << " |";
        std::cout << std::endl;
    }

private:
    std::string InsertThousand(size_t rValue, std::string rSeparator) const
    {
        std::string valueWithCommas = std::to_string(rValue);
        // length() is unsigned. avoid any subtraction with unsigned numbers ...
        int insertPosition = static_cast<int>(valueWithCommas.length()) - 3;
        while (insertPosition > 0)
        {
            valueWithCommas.insert(static_cast<unsigned>(insertPosition), rSeparator);
            insertPosition -= 3;
        }
        return valueWithCommas;
    }
};

class ReporterCSV : public Reporter
{
public:
    ReporterCSV(const std::string& mFileName)
    {
        mFile.open(mFileName);
    }

protected:
    virtual void WriteTestResult(const Result& rT) override
    {
        mFile << rT.mTestCase << "::" << rT.mTestName << ',' << rT.mTimeMean << std::endl;
    }

private:
    std::ofstream mFile;
};


class ReporterJUnit : public Reporter
{
public:
protected:
public:
    ReporterJUnit(const std::string& mFileName)
    {
        mFile.open(mFileName);
    }

protected:
    void WriteHeader() override
    {
        mFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        mFile << "<testsuites>\n";
        mFile << "  <testsuite name=\"All tests\">\n";
    }
    void WriteTestResult(const Result& rT) override
    {
        mFile << "    <testcase name=\"" << rT.mTestCase << "::" << rT.mTestName << "\" time=\"" << rT.mTimeMean
              << "\" \\>\n";
    }
    void WriteFooter() override
    {
        mFile << "  <\\testsuite>\n";
        mFile << "<\\testsuites>\n";
        mFile.flush();
    }

private:
    std::ofstream mFile;
};

//! @brief Simple Benchmark class that stores a name and a function to call
struct Benchmark
{
    std::string mTestCase;
    std::string mTestName;
    std::function<void(Runner&)> mFunction;
};

struct Warmup
{
    static void Test(BenchmarkInternal::Runner& runner)
    {
        while (runner.KeepRunningTime(0.5))
        {
        }
    }
};


class Holder
{
public:
    Holder()
        : mWarmup({"Warm", "up", &Warmup::Test})
    {
    }

    static Holder& GetInstance()
    {
        static Holder singleton;
        return singleton;
    }

    void Run(Reporter& rTestReporter)
    {
        Runner s;
        mWarmup.mFunction(s);
        for (auto& benchmark : mBenchmarks)
        {
            Runner s;
            benchmark.mFunction(s);
            rTestReporter.AddTestResult({benchmark.mTestCase, benchmark.mTestName, s.GetSeconds()});
        }
    }

    void AddBenchmark(Benchmark&& rBenchmark)
    {
        mBenchmarks.push_back(std::move(rBenchmark));
    }

private:
    Holder(const Holder&) = delete;
    Holder& operator=(const Holder&) = delete;
    std::vector<Benchmark> mBenchmarks;
    Benchmark mWarmup;
};


std::unique_ptr<Reporter> EvaluateArguments(int argc, char const* argv[])
{
    if (argc == 1) // default case
        return std::make_unique<ReporterCout>();

    std::string outputOption = argv[1];
    if (outputOption == "--cout")
    {
        return std::make_unique<ReporterCout>();
    }
    else if (outputOption == "--junit")
    {
        if (argc == 2)
        {
            std::cout << "You need to specify a file name! \n";
            return nullptr;
        }
        std::string outputFileName = argv[2];
        return std::make_unique<ReporterJUnit>(outputFileName);
    }
    else if (outputOption == "--csv")
    {
        if (argc == 2)
        {
            std::cout << "You need to specify a file name! \n";
            return nullptr;
        }
        std::string outputFileName = argv[2];
        return std::make_unique<ReporterCSV>(outputFileName);
    }
    std::cout << " Define the test output:\n";
    std::cout << " --cout              ... prints results to cout \n";
    std::cout << " --csv <FILE>.csv    ... prints results comma separated \n";
    std::cout << " --junit <FILE>.xml  ... prints results to FILE.xml \n";

    return nullptr;
}
} // namespace Benchmark

int main(int argc, char const* argv[])
{
    auto testReporter = BenchmarkInternal::EvaluateArguments(argc, argv);
    if (testReporter == nullptr)
        return -1;

    BenchmarkInternal::Holder::GetInstance().Run(*testReporter);
    testReporter->Report();
    return 0;
}

#define CONCAT(x, y) x##y
#define CLASS_NAME(x, y) CONCAT(x, y)
#define INSTANCE_NAME(x, line) CONCAT(x, line)


// ###################################################
//                BENCHMARK macros

#define BENCHMARK(TestCase, TestName, RunnerName)                                                                      \
    BENCHMARK_HELPER(CLASS_NAME(TestCase, TestName), INSTANCE_NAME(TestName, __LINE__), RunnerName, TestCase, TestName)

#define BENCHMARK_HELPER(NAME, INSTANCE, RUNNER, TestCase, TestName)                                                   \
    class NAME                                                                                                         \
    {                                                                                                                  \
    private:                                                                                                           \
        static void Test(BenchmarkInternal::Runner&);                                                                  \
                                                                                                                       \
    public:                                                                                                            \
        NAME()                                                                                                         \
        {                                                                                                              \
            BenchmarkInternal::Holder::GetInstance().AddBenchmark(                                                     \
                    BenchmarkInternal::Benchmark({#TestCase, #TestName, &Test}));                                      \
        }                                                                                                              \
    };                                                                                                                 \
    NAME INSTANCE;                                                                                                     \
    void NAME::Test(BenchmarkInternal::Runner& RUNNER)
// the actual function goes here (=just below BENCHMARK)
