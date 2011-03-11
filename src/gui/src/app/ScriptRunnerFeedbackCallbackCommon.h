/**\file
 * Common base class for ScriptRunner::FeedbackCallback implementations
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __SCRIPTRUNNERFEEDBACKCALLBACKCOMMON_H__
#define __SCRIPTRUNNERFEEDBACKCALLBACKCOMMON_H__

#include "ScriptRunner.h"

class ScriptRunnerFeedbackCallbackCommon : public nutogui::ScriptRunner::FeedbackCallback
{
public:
  void StartupComplete (bool success, const wxString& message) {}
  void Result (const nutogui::ResultDataSourceVTKPtr& result,
	       const wxString& caption) {}
  void ScriptOutput (nutogui::ScriptRunner::FeedbackCallback::OutputTarget target, const wxString& str) {}
  
  void ScriptRunStart () {}
  void ScriptRunEnd (bool success,
		     const nutogui::ScriptRunner::TracebackPtr& traceback) {}
};

#endif // __SCRIPTRUNNERFEEDBACKCALLBACKCOMMON_H__
