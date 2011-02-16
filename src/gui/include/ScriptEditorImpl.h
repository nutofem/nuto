#ifndef __NUTOGUI_SCRIPTEDITORIMPL_H__
#define __NUTOGUI_SCRIPTEDITORIMPL_H__

#include "GuiFrame.h"
#include "ScriptError.h"
#include "ScriptSource.h"
#include "TabCommonImpl.h"

#include <wx/wx.h>
#include <wx/docview.h>
#include <boost/weak_ptr.hpp>

class wxAuiToolBar;
class wxAuiToolBarEvent;
class wxStyledTextCtrl;
class wxStyledTextEvent;

#ifndef SCRIPTEDITOR_EXPORTED
  #define SCRIPTEDITOR_EXPORTED	NUTOGUI_IMPORT
#endif

namespace nutogui
{
  class SCRIPTEDITOR_EXPORTED ScriptEditorImpl : public TabCommonImpl,
						 public wxEvtHandler,
						 public ScriptSource,
						 public ScriptError
  {
    GuiFrame::TabCallbackWeakPtr callback;
    wxWindow* parentWindow;
    wxPanel* panel;
    wxStyledTextCtrl* editor;
    wxAuiToolBar* toolbar;
    wxStatusBar* statusBar;
    wxFont monospaceFont;
    wxStaticText* positionIndicator;
    
    wxMenu openPopup;
    wxFileHistory fileHistory;
    wxMenu savePopup;
    
    wxString currentHighlightWord;
    wxString newHighlightWord;
    void SetupEditorStyle ();
    void UpdateLineNumberMargin ();
    
    wxString fullPath;
    bool isChanged;
    wxString lastDir;
    
    int currentMarginDigits;
    int errorMarkerHandle;
    
    wxString GetDisplayFilename () const;
    void UpdateCaption ();
    bool DoFileOpen (const wxString& fullPath);
    bool DoFileSave (bool saveas);
    bool CanDiscardContents ();
    /// Save file history
    void SaveFilesConfig ();
    
    void OnFileNew (wxCommandEvent& event);
    void OnFileOpen (wxCommandEvent& event);
    void OnFileHistory (wxCommandEvent& event);
    void OnFileSave (wxCommandEvent& event);
    void OnFileSaveAs (wxCommandEvent& event);
    void OnToolbarDropdown (wxAuiToolBarEvent& event);
    void OnArtificialDropDown (wxCommandEvent& event);
    void DoDropDown (int toolID);
    
    void OnCopy (wxCommandEvent& event);
    void OnCut (wxCommandEvent& event);
    void OnPaste (wxCommandEvent& event);
    void OnUndo (wxCommandEvent& event);
    void OnRedo (wxCommandEvent& event);
    
    void OnUpdateUI (wxUpdateUIEvent& event);
    
    void OnEditorModified (wxStyledTextEvent& event);
    void OnEditorUpdateUI (wxStyledTextEvent& event);
    void UpdateStatusBar ();
    void OnEditorMarginClick (wxStyledTextEvent& event);
    void OnEditorCharAdded (wxStyledTextEvent& event);
    
    void OnSize (wxSizeEvent& event);
    void OnIdle (wxIdleEvent& event);
  public:
    ScriptEditorImpl ();
    ~ScriptEditorImpl ();
    
    wxWindow* CreateContents (const GuiFrame::TabCallbackWeakPtr& callback,
			      wxWindow* parentWindow);
			      
    void TabActivate ();
    
    /**\name ScriptSource implementation
     * @{ */
    wxString GetSourceName() const;
    wxString GetSourceString() const;
    /** @} */
    
    /**\name ScriptError implementation
     * @{ */
    void ClearErrorLocation ();
    bool SetErrorLocation (const ScriptRunner::TracebackPtr& traceback);
    void DisplayErrorLocation ();
    /** @} */
			      
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_SCRIPTEDITORIMPL_H__
