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

#include "ResultDataSourceVTKFromFile.h"
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

void ScriptResultDispatcher::Result (const wxString& resultTempFile,
				     const wxString& caption,
				     const wxString& resultName)
{
  ResultDataSourceVTKFromFilePtr resultsContainer;
  {
    ResultsMap::const_iterator existingResult = results.find (caption);
    if (existingResult != results.end())
      resultsContainer = existingResult->second;
    else
    {
      resultsContainer = boost::make_shared<ResultDataSourceVTKFromFile> ();
      results[caption] = resultsContainer;
    }
  }
  resultsContainer->AddDataSet (resultTempFile.fn_str(), resultName);
  wxRemoveFile (resultTempFile);
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
  
  // Create result tabs for collected results
  for (ResultsMap::const_iterator result = results.begin();
       result != results.end();
       ++result)
  {
    // Create a result frame
    wxString captionStr (result->first);
    if (captionStr.IsEmpty())
      captionStr = wxT ("Results");
    
    wxString fullCaption (wxString::Format (wxT ("%s (%s)"),
					    captionStr.c_str(),
					    scriptSource->GetSourceName().c_str()));
    
    nutogui::GuiFrame::TabPtr resultViewerTab =
      boost::make_shared<nutogui::ResultViewerImpl> (result->second, fullCaption);
    frame->AddTab (resultViewerTab, true);
  }
  // Clean up all collected results
  results.clear();
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
