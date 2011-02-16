#ifndef __APP_SCRIPTCONTROLBAR_H__
#define __APP_SCRIPTCONTROLBAR_H__

#include "GuiFrame.h"
#include "ScriptResultDispatcher.h"
#include "ScriptRunner.h"
#include "ScriptRunnerFeedbackCallbackCommon.h"
#include "ScriptSource.h"

#include <wx/wx.h>
#include <boost/enable_shared_from_this.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>

class wxAuiToolBar;
class wxAuiToolBarItem;

class ScriptControlBar : public nutogui::GuiFrame::Pane,
			 public wxEvtHandler
{
  nutogui::ScriptRunnerPtr scriptRunner;
  nutogui::ScriptSourcePtr scriptSource;
  ScriptResultDispatcherPtr resultDispatch;
  wxWindow* parentWindow;
  wxAuiToolBar* bar;
  wxAuiToolBarItem* runningIndicator;
  wxString currentScriptLine;
  enum { Normal = 0, Error = 1 };
  int scriptOutputLevel;
  bool scriptRunnerUseable;

  // Separate class to avoid circular ref
  class FeedbackCallbackImpl : public ScriptRunnerFeedbackCallbackCommon
  {
    ScriptControlBar* owner;
  public:
    FeedbackCallbackImpl (ScriptControlBar* owner)
     : owner (owner) {}
    
    /**\name nutogui::ScriptRunner::StartupCallback implementation
    * @{ */
    void StartupComplete (bool success, const wxString& message)
    {
      owner->StartupComplete (success, message);
    }
    void ScriptOutput (nutogui::ScriptRunner::FeedbackCallback::OutputTarget target,
		       const wxString& str)
    {
      owner->ScriptOutput (target, str);
    }
    void ScriptRunStart ()
    {
      owner->ScriptRunStart();
    }
    void ScriptRunEnd (bool success,
		       const nutogui::ScriptRunner::TracebackPtr& traceback)
    {
      owner->ScriptRunEnd (success);
    }
    /** @} */
  };
  
  void StartupComplete (bool success, const wxString& message);
  void ScriptOutput (nutogui::ScriptRunner::FeedbackCallback::OutputTarget target, const wxString& str);
  void ScriptRunStart ();
  void ScriptRunEnd (bool success);
  void OnRunScript (wxCommandEvent& event);
  void OnGoToError (wxCommandEvent& event);
  void OnGoToErrorUIUpdate (wxUpdateUIEvent& event);
  void OnShowOutput (wxCommandEvent& event);
  
  void EnsureRunningIndicator ();
  void SetIndicatorLabel (const wxString& label);
  void SetIndicatorTooltip (const wxString& tip);
  std::vector<wxBitmap> throbberFrames;
  bool throbbing;
  wxTimer throbberTimer;
  size_t throbberFrame;
  void LoadThrobber ();
  void StartThrobber ();
  void StopThrobber ();
  void OnThrobberTimer (wxTimerEvent&);
public:
  ScriptControlBar (const nutogui::ScriptRunnerPtr& scriptRunner,
		    const nutogui::ScriptSourcePtr& scriptSource,
		    const ScriptResultDispatcherPtr& resultDispatch);
  
  /**\name nutogui::GuiFrame::Pane implementation
   * @{ */
  wxWindow* CreateContents (const nutogui::GuiFrame::PaneCallbackWeakPtr& callback,
			    wxWindow* parentWindow);
			    
  bool IsToolbar() const { return true; }
  /** @} */
  
  DECLARE_EVENT_TABLE()
};

#endif // __APP_SCRIPTCONTROLBAR_H__
