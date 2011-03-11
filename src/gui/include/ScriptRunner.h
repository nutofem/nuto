/**\file
 * Interface for script runner.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_SCRIPTRUNNER_H__
#define __NUTOGUI_SCRIPTRUNNER_H__

#include "ResultDataSourceVTK.h"

#include <wx/string.h>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace nutogui
{
  struct ScriptRunner
  {
    /// Entry in a traceback
    struct TracebackEntry
    {
      wxString module;
      wxString filename;
      int line;
    };
    typedef std::vector<TracebackEntry> Traceback;
    typedef boost::shared_ptr<const Traceback> TracebackPtr;
    
    /// Interface to receive script state notifications.
    struct FeedbackCallback
    {
      virtual ~FeedbackCallback() {}
      
      virtual void StartupComplete (bool success, const wxString& message) = 0;
      virtual void Result (const wxString& resultTempFile,
			   const wxString& caption,
			   const wxString& resultName) = 0;

      enum OutputTarget { StdOut, StdErr };
      virtual void ScriptOutput (OutputTarget target, const wxString& str) = 0;
      
      virtual void ScriptRunStart () = 0;
      virtual void ScriptRunEnd (bool success, const TracebackPtr& traceback) = 0;
    };
    typedef boost::shared_ptr<FeedbackCallback> FeedbackCallbackPtr;
    
    virtual ~ScriptRunner() {}
    
    virtual void AddFeedbackCallback (const FeedbackCallbackPtr& callback) = 0;
    
    virtual void StartScript (const wxString& source) = 0;
  };
  typedef boost::shared_ptr<ScriptRunner> ScriptRunnerPtr;
} // namespace nutogui

#endif // __NUTOGUI_SCRIPTRUNNER_H__
