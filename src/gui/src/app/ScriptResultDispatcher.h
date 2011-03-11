/**\file
 * Show results from a script run
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __SCRIPTRESULTDISPATCHER_H__
#define __SCRIPTRESULTDISPATCHER_H__

#include "GuiFrame.h"
#include "ScriptError.h"
#include "ScriptSource.h"

#include "ScriptRunnerFeedbackCallbackCommon.h"

class ScriptOutputTab;

class ScriptResultDispatcher : public ScriptRunnerFeedbackCallbackCommon
{
  nutogui::GuiFrame* frame;
  nutogui::ScriptSourcePtr scriptSource;
  nutogui::ScriptErrorPtr scriptError;
  boost::shared_ptr<ScriptOutputTab> scriptOutput;
  
  bool scriptHadWarnings;
  bool hasErrorLocation;
public:
  ScriptResultDispatcher (nutogui::GuiFrame* frame,
			  const nutogui::ScriptSourcePtr& scriptSource,
			  const nutogui::ScriptErrorPtr& scriptError)
   : frame (frame), scriptSource (scriptSource), scriptError (scriptError),
     hasErrorLocation (false)
  {}
  
  /**\name nutogui::ScriptRunner::FeedbackCallback implementation
   * @{ */
  void StartupComplete (bool success, const wxString& message);
  void Result (const nutogui::ResultDataSourceVTKPtr& result,
	       const wxString& caption);
		       
  void ScriptOutput (OutputTarget target, const wxString& str);
  
  void ScriptRunStart ();
  void ScriptRunEnd (bool success,
		     const nutogui::ScriptRunner::TracebackPtr& traceback);
  /** @} */
  
  void ShowScriptOutput ();
  bool CanGoToErrorLocation ();
  void GoToErrorLocation ();
};
typedef boost::shared_ptr<ScriptResultDispatcher> ScriptResultDispatcherPtr;

#endif // __SCRIPTRESULTDISPATCHER_H__
