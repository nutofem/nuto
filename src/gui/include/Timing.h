/**\file
 * Timing helper
 */
/*
 * Written 2009 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __TIMING_H__
#define __TIMING_H__

#include <wx/stopwatch.h>
#include <wx/log.h>

namespace nutogui
{
  /**
   * Helper class to time how long a block takes to process.
   * The time between construction and destruction of the timing instance
   * is measured and reported on destruction.
   */
  class Timing
  {
    wxStopWatch sw;
    typedef long TimeType;
#define __TIMING_FORMAT		"%ld ms"
    TimeType GetTime() const { return sw.Time(); }

    const char* description;
    TimeType startTime;
  public:
    Timing (const char* descr) : description (descr)
    {
      startTime = GetTime();
    }
    ~Timing()
    {
      TimeType endTime = GetTime();
      wxLogDebug (wxT ("%s: " __TIMING_FORMAT), description, endTime-startTime);
    }
    /// Print out the time that passed so far
    void PrintIntermediate (const char* descr) const
    {
      TimeType time = GetTime();
      wxLogDebug (wxT ("%s: %s " __TIMING_FORMAT), description, descr, time-startTime);
    }
#undef __TIMING_FORMAT
  };
} // namespace spectacles

#endif // __TIMING_H__
