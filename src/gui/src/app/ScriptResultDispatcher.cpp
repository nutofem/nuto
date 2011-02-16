/**\file
 * Show results from a script run
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "ScriptResultDispatcher.h"

#include "uicommon/Quote.h"
#include "ResultViewerImpl.h"

#include "ScriptOutputTab.h"

#include <boost/make_shared.hpp>
#include <wx/wx.h>
#include <wx/tokenzr.h>
#include <boost/concept_check.hpp>

void ScriptResultDispatcher::StartupComplete (bool success, const wxString& message)
{
  if (success)
  {
    scriptOutput = boost::make_shared<ScriptOutputTab> ();
    frame->AddTab (scriptOutput, false);
  }
}

void ScriptResultDispatcher::ResultDataFile (const wxString& fileName,
					     const wxString& caption)
{
  // Create a result frame
  wxString captionStr (caption);
  if (captionStr.IsEmpty())
    captionStr = wxT ("Results");
  
  wxString fullCaption (wxString::Format (wxT ("%s (%s)"),
					  captionStr.c_str(),
					  scriptSource->GetSourceName().c_str()));
  
  nutogui::GuiFrame::TabPtr resultViewerTab =
    boost::make_shared<nutogui::ResultViewerImpl> (fileName, fullCaption, true);
  frame->AddTab (resultViewerTab, true);
}

void ScriptResultDispatcher::ScriptOutput (OutputTarget target, const wxString& str)
{
  // Log output: one log entry for each line, strip newlines
  wxStringTokenizer tokenize (str, wxT("\n"));
  wxString line;
  while (tokenize.HasMoreTokens())
  {
    line = tokenize.GetNextToken();
    if (target == StdOut)
    {
      wxLogMessage (wxT ("Python script output: %s"), line.c_str());
    }
    else
    {
      wxLogWarning (wxT ("Python script error: %s"), line.c_str());
    }
  }
  // Script output tab: output verbatim
  if (target == StdOut)
  {
    scriptOutput->AppendOutput (str);
  }
  else
  {
    scriptOutput->AppendErrorOutput (str);
    scriptHadWarnings = true;
  }
}

void ScriptResultDispatcher::ScriptRunStart ()
{
  scriptHadWarnings = false;
  scriptOutput->AppendNote (wxString::Format (wxT ("Script %s started: %s\n"),
					      uicommon::Quote::Single (scriptSource->GetSourceName()).c_str(),
					      wxDateTime::Now().Format().c_str()));
}

void ScriptResultDispatcher::ScriptRunEnd (bool success,
					   const nutogui::ScriptRunner::TracebackPtr& traceback)
{
  wxString scriptState;
  if (!success)
    scriptState = wxT ("with exception");
  else if (scriptHadWarnings)
    scriptState = wxT ("with warnings or errors");
  else
    scriptState = wxT ("successful");
  scriptOutput->AppendNote (wxString::Format (wxT ("Script %s finished: %s (%s)\n"),
					      uicommon::Quote::Single (scriptSource->GetSourceName()).c_str(),
					      wxDateTime::Now().Format().c_str(),
					      scriptState.c_str()));
					      
  if (success)
  {
    scriptError->ClearErrorLocation ();
    hasErrorLocation = false;
  }
  else
  {
    hasErrorLocation = scriptError->SetErrorLocation (traceback);
  }
}
  
void ScriptResultDispatcher::ShowScriptOutput ()
{
  frame->ActivateTab (scriptOutput.get());
}

bool ScriptResultDispatcher::CanGoToErrorLocation ()
{
  return hasErrorLocation;
}

void ScriptResultDispatcher::GoToErrorLocation ()
{
  scriptError->DisplayErrorLocation ();
}
