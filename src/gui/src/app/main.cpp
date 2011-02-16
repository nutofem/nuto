#include "common.h"

#include "main.h"

#include "uicommon/ArtProvider.h"
#include "GuiFrameImpl.h"
#include "ScriptControlBar.h"
#include "ScriptEditorImpl.h"
#include "ScriptResultDispatcher.h"

#include "ScriptRunner/ScriptRunnerImpl.h"

#include <boost/make_shared.hpp>

IMPLEMENT_APP(Application)

static const wxChar appName[] = wxT ("NuToGUI");

bool Application::OnInit()
{
#if wxABI_VERSION >= 20900
  SetAppDisplayName (appName);
#endif
  SetAppName (appName);
  SetVendorName (wxT ("BU Weimar ISM"));

  if (!wxApp::OnInit())
    return false;
  
  // Set up logging targets
  wxLogStream* log_cerr = new wxLogStream;
  wxLog::SetActiveTarget (log_cerr);
  
  // Set up art provider (for custom + fallback art)
  wxArtProvider::Push (new uicommon::ArtProvider);

  nutogui::GuiFrameImpl* frame = new nutogui::GuiFrameImpl (appName);
  
  wxLogWindow* log_window = new wxLogWindow (frame, wxT ("Message Log"), false);

  nutogui::ScriptSourcePtr scriptSource;
  nutogui::ScriptErrorPtr scriptError;
  // Show script editor by default
  {
    boost::shared_ptr<nutogui::ScriptEditorImpl> scriptEditorTab = boost::make_shared<nutogui::ScriptEditorImpl> ();
    frame->AddTab (scriptEditorTab);
    scriptSource = scriptEditorTab;
    scriptError = scriptEditorTab;
  }
  // Set up script runner
  nutogui::ScriptRunnerPtr scriptRunner = boost::make_shared<ScriptRunnerImpl> ();
  // Class to display script results, when they arrive
  ScriptResultDispatcherPtr resultDispatch =
    boost::make_shared<ScriptResultDispatcher> (frame, scriptSource, scriptError);
  scriptRunner->AddFeedbackCallback (resultDispatch);
  // Toolbar to launch scripts
  {
    nutogui::GuiFrame::PanePtr scriptControlPane =
      boost::make_shared<ScriptControlBar> (scriptRunner, scriptSource, resultDispatch);
    frame->AddPane (wxT ("scriptControlBar"), nutogui::GuiFrame::Bottom,
		    scriptControlPane);
  }

  frame->SetLogWindow (log_window);
  frame->Show (true);
  SetTopWindow (frame);
  return true;
}
