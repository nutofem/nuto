#include "common.h"

#include "ScriptControlBar.h"

#include "uicommon/ArtProvider.h"
#include "uicommon/MessageDialog.h"

#include <wx/wx.h>
#include <wx/artprov.h>
#include <wx/aui/auibar.h>
#include <wx/aui/framemanager.h>
#include <wx/mstream.h>

#include <boost/make_shared.hpp>

enum
{
  ID_RunScript = 100,
  ID_GoToError,
  ID_ShowOutput,
  
  TIMER_Throbber
};

BEGIN_EVENT_TABLE(ScriptControlBar, wxEvtHandler)
  EVT_TIMER(TIMER_Throbber, ScriptControlBar::OnThrobberTimer)
END_EVENT_TABLE()

ScriptControlBar::ScriptControlBar (const nutogui::ScriptRunnerPtr& scriptRunner,
				    const nutogui::ScriptSourcePtr& scriptSource,
				    const ScriptResultDispatcherPtr& resultDispatch)
 : scriptRunner (scriptRunner), scriptSource (scriptSource), resultDispatch (resultDispatch),
   parentWindow (parentWindow), bar (nullptr), runningIndicator (nullptr),
   scriptRunnerUseable (false),
   throbbing (false), throbberTimer (this, TIMER_Throbber)
{
  assert (scriptRunner);
  assert (scriptSource);
  
  nutogui::ScriptRunner::FeedbackCallbackPtr callback (
    boost::make_shared<FeedbackCallbackImpl> (this));
  scriptRunner->AddFeedbackCallback (callback);
}

wxWindow* ScriptControlBar::CreateContents (const nutogui::GuiFrame::PaneCallbackWeakPtr& callback_,
					    wxWindow* parentWindow)
{
  nutogui::GuiFrame::PaneCallbackPtr callback (callback_);
  callback->SetCaption (wxT ("Script Control"));
  callback->AddAccelerator (wxAcceleratorEntry (0, WXK_F9, ID_RunScript));
  this->parentWindow = parentWindow;
  
  bar = new wxAuiToolBar (parentWindow, wxID_ANY,
			  wxDefaultPosition, wxDefaultSize,
			  wxAUI_TB_HORZ_TEXT);
  bar->SetToolBitmapSize (wxArtProvider::GetSizeHint (wxART_TOOLBAR));
  bar->AddTool (ID_RunScript, wxT ("Run script"),
		wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(system-run), wxART_TOOLBAR),
		wxT ("Run the current script file (F9)"));
  bar->EnableTool (ID_RunScript, false);
  bar->Realize();
  
  /* We have to explicitly connect to events instead of using an event table
     as an event handler pushed onto the bar doesn't actually get events from
     it :/ */
  bar->Connect (ID_RunScript, wxEVT_COMMAND_MENU_SELECTED,
		wxCommandEventHandler (ScriptControlBar::OnRunScript),
		nullptr, this);
  bar->Connect (ID_GoToError, wxEVT_COMMAND_MENU_SELECTED,
		wxCommandEventHandler (ScriptControlBar::OnGoToError),
		nullptr, this);
  bar->Connect (ID_GoToError, wxEVT_UPDATE_UI,
		wxUpdateUIEventHandler (ScriptControlBar::OnGoToErrorUIUpdate),
		nullptr, this);
  bar->Connect (ID_ShowOutput, wxEVT_COMMAND_MENU_SELECTED,
		wxCommandEventHandler (ScriptControlBar::OnShowOutput),
		nullptr, this);
		
  LoadThrobber();
  
  return bar;
}

void ScriptControlBar::StartupComplete (bool success, const wxString& message)
{
  scriptRunnerUseable = success;
  if (success)
  {
    bar->EnableTool (ID_RunScript, true);
    bar->Realize();
  }
  else
  {
    uicommon::MessageDialog dlg;
    dlg.SetIcon (wxICON_ERROR);
    dlg.SetMessageHeading (wxT ("Could not start up Python."));
    dlg.SetMessageBody (wxString::Format (wxT ("Running scripts is not possible.\n\nPython output was:\n%s"),
					  message.c_str()));
    dlg.OK (parentWindow);
  }
}

void ScriptControlBar::ScriptOutput (nutogui::ScriptRunner::FeedbackCallback::OutputTarget target,
				     const wxString& str)
{
  if (target == nutogui::ScriptRunner::FeedbackCallback::StdErr)
    scriptOutputLevel = Error;
  
  currentScriptLine.Append (str);
  int lfPos = currentScriptLine.Find ('\n');
  wxString displayLine;
  while (lfPos != wxNOT_FOUND)
  {
    wxString newPart = currentScriptLine.Mid (lfPos+1);
    if (newPart.IsEmpty())
    {
      displayLine = currentScriptLine.Mid (0, lfPos);
      // Stash a newline in currentScriptLine so the next output will provoke a change of display
      currentScriptLine = wxT("\n");
      break;
    }
    else
    {
      displayLine = newPart;
      currentScriptLine = newPart;
      lfPos = currentScriptLine.Find ('\n');
    }
  }
  
  if (!displayLine.IsEmpty())
  {
    SetIndicatorLabel (displayLine);
  }
}

void ScriptControlBar::ScriptRunStart ()
{
  scriptOutputLevel = Normal;
  currentScriptLine.Empty();
  
  EnsureRunningIndicator();
  StartThrobber();
  SetIndicatorTooltip (wxT ("Script is running..."));
  SetIndicatorLabel (wxEmptyString);
  
  bar->EnableTool (ID_RunScript, false);
  bar->Realize();
}

void ScriptControlBar::ScriptRunEnd (bool success)
{
  StopThrobber ();
  
  wxBitmap indicatorBitmap;
  wxString indicatorMsg;
  if (!success)
  {
    indicatorBitmap = wxArtProvider::GetBitmap (wxART_ERROR, wxART_TOOLBAR);
    indicatorMsg = wxT ("Script failed with an error.");
  }
  else if (scriptOutputLevel == Error)
  {
    indicatorBitmap = wxArtProvider::GetBitmap (wxART_WARNING, wxART_TOOLBAR);
    indicatorMsg = wxT ("Script ran successfully, but there are warning or error messages.");
  }
  else
  {
    indicatorBitmap = wxArtProvider::GetBitmap (wxART_INFORMATION, wxART_TOOLBAR);
    indicatorMsg = wxT ("Script ran successfully.");
  }
  runningIndicator->SetBitmap (indicatorBitmap);
  SetIndicatorTooltip (indicatorMsg);
  
  bar->EnableTool (ID_RunScript, true);
  bar->Realize();
}

void ScriptControlBar::OnRunScript (wxCommandEvent& event)
{
  if (!scriptRunnerUseable) return;
  
  scriptRunner->StartScript (scriptSource->GetSourceString());
}

void ScriptControlBar::OnGoToError (wxCommandEvent& event)
{
  resultDispatch->GoToErrorLocation ();
}

void ScriptControlBar::OnGoToErrorUIUpdate (wxUpdateUIEvent& event)
{
  event.Enable (resultDispatch->CanGoToErrorLocation ());
}

void ScriptControlBar::OnShowOutput (wxCommandEvent& event)
{
  resultDispatch->ShowScriptOutput();
}

void ScriptControlBar::EnsureRunningIndicator ()
{
  if (runningIndicator) return;
  
  bar->AddSeparator();
  bar->AddTool (ID_GoToError, wxEmptyString,
		wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(goto-error), wxART_TOOLBAR),
		wxT ("Go to error location"));
  bar->EnableTool (ID_GoToError, false);
  bar->AddTool (ID_ShowOutput, wxEmptyString,
		throbberFrames[1]);
  runningIndicator = bar->FindToolByIndex (bar->GetToolCount() - 1);
}

void ScriptControlBar::SetIndicatorLabel (const wxString& label)
{
  runningIndicator->SetLabel (label);
  bar->Realize();
  
  // Toolbar size might have changed, so tell manager about it
  wxAuiManager* auimgr = wxAuiManager::GetManager (bar);
  if (auimgr)
  {
    wxAuiPaneInfo& paneInfo = auimgr->GetPane (bar);
    paneInfo.BestSize (bar->GetSize ());
    auimgr->Update();
  }
}

void ScriptControlBar::SetIndicatorTooltip (const wxString& tip)
{
  runningIndicator->SetShortHelp (wxString::Format (wxT ("%s Click to see script output."),
						    tip.c_str()));
}

#include "throbber.png.inc"

void ScriptControlBar::LoadThrobber ()
{
  wxSize throbberFrameSize (22, 22);
  wxSize throbberTargetSize (wxArtProvider::GetSizeHint (wxART_TOOLBAR));
  
  wxMemoryInputStream throbberStream (throbber_png_data, throbber_png_size);
  wxImage throbberStrip (throbberStream, wxT("image/png"));
  
  int rows = throbberStrip.GetHeight() / throbberFrameSize.GetHeight();
  int cols = throbberStrip.GetWidth() / throbberFrameSize.GetWidth();
  
  int curRow = 0, curCol = 0;
  for (int n = 0; n < (rows*cols); n++)
  {
    /* Extract a single frame from the strip */
    wxRect subRect (curCol * throbberFrameSize.GetWidth(),
		    curRow * throbberFrameSize.GetHeight(),
		    throbberFrameSize.GetWidth(),
		    throbberFrameSize.GetHeight());
    wxImage throbberFrame (throbberStrip.GetSubImage (subRect));
    curCol += 1;
    if (curCol >= cols)
    {
      curCol = 0;
      curRow += 1;
    }

    // Resize frame to toolbar size
    // @@@ FIXME: downsample when TB size is smaller
    throbberFrame.Resize (throbberTargetSize,
			  wxPoint ((throbberTargetSize.GetWidth() - throbberFrameSize.GetWidth()) / 2,
				   (throbberTargetSize.GetHeight() - throbberFrameSize.GetHeight()) / 2));
				   
    throbberFrames.push_back (wxBitmap (throbberFrame));
  }
}

void ScriptControlBar::StartThrobber ()
{
  runningIndicator->SetBitmap (throbberFrames[1]);
  throbberFrame = 0;
  throbbing = true;
  
  const int throbberCycleTime = 2000;
  throbberTimer.Start (throbberCycleTime / (throbberFrames.size()-1));
}

void ScriptControlBar::StopThrobber ()
{
  throbbing = false;
  throbberTimer.Stop ();
}

void ScriptControlBar::OnThrobberTimer (wxTimerEvent&)
{
  if (!throbbing) return;
  
  throbberFrame = (throbberFrame+1)%(throbberFrames.size()-1);
  runningIndicator->SetBitmap (throbberFrames[throbberFrame+1]);
  bar->Refresh();
}
