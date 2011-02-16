#include "common.h"
#include "ScriptEditorImpl.h"

#include "uicommon/MessageDialog.h"
#include "uicommon/Quote.h"

#include <wx/wx.h>
#include <wx/artprov.h>
#include <wx/aui/auibar.h>
#include <wx/config.h>
#include <wx/file.h>
#include <wx/filename.h>
#include <wx/stc/stc.h>
#include <wx/stockitem.h>

#include <boost/foreach.hpp>
#include <vector>

namespace nutogui
{
  using uicommon::Quote;
  
  enum
  {
    marginLineNumber = 0,
    marginFolding = 2
  };
  
  enum
  {
    ID_Editor = 1,
    ID_Toolbar,
    
    ID_DropDownOpen
  };
  
  static const int defaultLineNumberDigits = 3;
  
  BEGIN_EVENT_TABLE(ScriptEditorImpl, wxEvtHandler)
    EVT_MENU(wxID_NEW, ScriptEditorImpl::OnFileNew)
    EVT_MENU(wxID_OPEN, ScriptEditorImpl::OnFileOpen)
    EVT_AUITOOLBAR_TOOL_DROPDOWN(wxID_OPEN, ScriptEditorImpl::OnToolbarDropdown)
    EVT_MENU_RANGE(wxID_FILE1, wxID_FILE1+9, ScriptEditorImpl::OnFileHistory)
    EVT_MENU(wxID_SAVE, ScriptEditorImpl::OnFileSave)
    EVT_MENU(wxID_SAVEAS, ScriptEditorImpl::OnFileSaveAs)
    EVT_AUITOOLBAR_TOOL_DROPDOWN(wxID_SAVE, ScriptEditorImpl::OnToolbarDropdown)
    EVT_MENU(ID_DropDownOpen, ScriptEditorImpl::OnArtificialDropDown)
    
    EVT_MENU(wxID_COPY, ScriptEditorImpl::OnCopy)
    EVT_UPDATE_UI(wxID_COPY, ScriptEditorImpl::OnUpdateUI)
    EVT_MENU(wxID_CUT, ScriptEditorImpl::OnCut)
    EVT_UPDATE_UI(wxID_CUT, ScriptEditorImpl::OnUpdateUI)
    EVT_MENU(wxID_PASTE, ScriptEditorImpl::OnPaste)
    EVT_UPDATE_UI(wxID_PASTE, ScriptEditorImpl::OnUpdateUI)
    EVT_MENU(wxID_UNDO, ScriptEditorImpl::OnUndo)
    EVT_UPDATE_UI(wxID_UNDO, ScriptEditorImpl::OnUpdateUI)
    EVT_MENU(wxID_REDO, ScriptEditorImpl::OnRedo)
    EVT_UPDATE_UI(wxID_REDO, ScriptEditorImpl::OnUpdateUI)
    
    EVT_STC_MODIFIED(ID_Editor, ScriptEditorImpl::OnEditorModified)
    EVT_STC_UPDATEUI(ID_Editor, ScriptEditorImpl::OnEditorUpdateUI)
    EVT_STC_MARGINCLICK(ID_Editor, ScriptEditorImpl::OnEditorMarginClick)
    EVT_STC_CHARADDED(ID_Editor, ScriptEditorImpl::OnEditorCharAdded)
    
    EVT_SIZE(ScriptEditorImpl::OnSize)
    EVT_IDLE(ScriptEditorImpl::OnIdle)
  END_EVENT_TABLE()

  ScriptEditorImpl::ScriptEditorImpl ()
   : parentWindow (nullptr), panel (nullptr), isChanged (false),
     currentMarginDigits (0),
     errorMarkerHandle (-1)
  {}

  ScriptEditorImpl::~ScriptEditorImpl()
  {
    if (panel) panel->PopEventHandler ();
  }
  
  typedef std::vector<wxAcceleratorEntry> AcceleratorVector;
  
  static void AddStockAccel (AcceleratorVector& accelEntries, int id)
  {
    wxAcceleratorEntry accel (wxGetStockAccelerator (id));
    if (accel.IsOk())
      accelEntries.push_back (accel);
  }

  wxWindow* ScriptEditorImpl::CreateContents (const GuiFrame::TabCallbackWeakPtr& callback,
					      wxWindow* parentWindow)
  {
    this->callback = callback;
    this->parentWindow = parentWindow;
    
    UpdateCaption ();
    
    wxFont defFont (wxSystemSettings::GetFont (wxSYS_DEFAULT_GUI_FONT));
    monospaceFont = wxFont (defFont.GetPointSize(), wxFONTFAMILY_MODERN, defFont.GetStyle(), defFont.GetWeight());
    
    panel = new wxPanel (parentWindow);
    panel->PushEventHandler (this);
    wxSizer* sizer = new wxBoxSizer (wxVERTICAL);
    
    wxString openDropDownAccelStr (wxT ("Alt+F"));
    
    toolbar = new wxAuiToolBar (panel, ID_Toolbar, wxDefaultPosition, wxDefaultSize,
				wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
    toolbar->SetToolBitmapSize (wxArtProvider::GetSizeHint (wxART_TOOLBAR));
    toolbar->AddTool (wxID_OPEN, wxT ("Open"),
		      wxArtProvider::GetBitmap (wxART_FILE_OPEN, wxART_TOOLBAR),
		      wxString::Format (wxT ("Open a script file (%s)"),
					openDropDownAccelStr.c_str()));
    toolbar->SetToolDropDown (wxID_OPEN, true); // Drop down: file history
    toolbar->AddTool (wxID_SAVE, wxT ("Save"),
		      wxArtProvider::GetBitmap (wxART_FILE_SAVE, wxART_TOOLBAR),
		      wxT ("Save script file"));
    toolbar->SetToolDropDown (wxID_SAVE, true); // Drop down: menu with "Save as"
    toolbar->AddSeparator();
    toolbar->AddTool (wxID_COPY, wxT ("Copy"),
		      wxArtProvider::GetBitmap (wxART_COPY, wxART_TOOLBAR),
		      wxT ("Copy to clipboard"));
    toolbar->AddTool (wxID_CUT, wxT ("Cut"),
		      wxArtProvider::GetBitmap (wxART_CUT, wxART_TOOLBAR),
		      wxT ("Cut to clipboard"));
    toolbar->AddTool (wxID_PASTE, wxT ("Paste"),
		      wxArtProvider::GetBitmap (wxART_PASTE, wxART_TOOLBAR),
		      wxT ("Paste from clipboard"));
    toolbar->AddSeparator();
    toolbar->AddTool (wxID_UNDO, wxT ("Undo"),
		      wxArtProvider::GetBitmap (wxART_UNDO, wxART_TOOLBAR),
		      wxT ("Undo last change"));
    toolbar->AddTool (wxID_REDO, wxT ("Redo"),
		      wxArtProvider::GetBitmap (wxART_REDO, wxART_TOOLBAR),
		      wxT ("Redo undone change"));
    toolbar->Realize();
    sizer->Add (toolbar, 0, wxEXPAND);
    
    editor = new wxStyledTextCtrl (panel, ID_Editor);
    SetupEditorStyle ();
    UpdateLineNumberMargin ();
    sizer->Add (editor, 1, wxEXPAND);
    
    statusBar = new wxStatusBar (panel, wxID_ANY, 0);
    statusBar->SetFieldsCount (1);
    static const int styles[] = { wxSB_FLAT };
    statusBar->SetStatusStyles (1, styles);
    sizer->Add (statusBar, 0, wxEXPAND);
    
    positionIndicator = new wxStaticText (statusBar, wxID_ANY, wxEmptyString,
					  wxDefaultPosition, wxDefaultSize);
    positionIndicator->SetFont (monospaceFont);
    
    panel->SetSizer (sizer);
    
    // Keyboard shortcuts
    std::vector<wxAcceleratorEntry> accelEntries;
    
    {
      openPopup.Append (wxID_NEW);
      AddStockAccel (accelEntries, wxID_NEW);
    }
    {
      wxMenuItem* openItem = new wxMenuItem(&openPopup, wxID_OPEN,
					    wxT ("&Open...\tCtrl+O"));
      openItem->SetBitmap (wxArtProvider::GetBitmap (wxART_FILE_OPEN, wxART_MENU));
      openPopup.Append (openItem);
      AddStockAccel (accelEntries, wxID_OPEN);
    }
    fileHistory.UseMenu (&openPopup);
    
    {
      wxMenuItem* saveItem = new wxMenuItem(&savePopup, wxID_SAVE);
      saveItem->SetBitmap (wxArtProvider::GetBitmap (wxART_FILE_SAVE, wxART_MENU));
      savePopup.Append (saveItem);
      AddStockAccel (accelEntries, wxID_SAVE);
    }
    {
      wxMenuItem* saveAsItem = new wxMenuItem(&savePopup, wxID_SAVEAS,
					      wxT ("Save &As...\tShift+Ctrl+S"));
      savePopup.Append (saveAsItem);
      wxAcceleratorEntry* accel = saveAsItem->GetAccel();
      if (accel && accel->IsOk())
	accelEntries.push_back (wxAcceleratorEntry (accel->GetFlags(), accel->GetKeyCode(), wxID_SAVEAS));
      delete accel;
    }
    
    // Add hot key to pop up open drop down. Mostly useful for accessing history items
    {
      wxAcceleratorEntry openDropDownAccel (0, 0, ID_DropDownOpen);
      if (openDropDownAccel.FromString (wxString (wxT ("\t")) + openDropDownAccelStr))
      {
	accelEntries.push_back (openDropDownAccel);
      }
    }
    
    wxAcceleratorTable accelTable (accelEntries.size(), accelEntries.data());
    panel->SetAcceleratorTable (accelTable);
    
    {
      wxConfigBase& config = *(wxConfigBase::Get ());
      config.SetPath (wxT ("ScriptEditor"));
      fileHistory.Load (config);
      config.Read (wxT("lastDirectory"), &lastDir);
      config.SetPath (wxT (".."));
    }
    
    return panel;
  }

  void ScriptEditorImpl::TabActivate ()
  {
    editor->SetFocus ();
  }

  wxString ScriptEditorImpl::GetSourceName() const
  {
    return GetDisplayFilename ();
  }

  wxString ScriptEditorImpl::GetSourceString() const
  {
    return editor->GetText();
  }
  
  enum
  {
    MARKER_Error = 0
  };

  void ScriptEditorImpl::ClearErrorLocation ()
  {
    editor->MarkerDeleteHandle (errorMarkerHandle);
    errorMarkerHandle = -1;
  }

  bool ScriptEditorImpl::SetErrorLocation (const ScriptRunner::TracebackPtr& traceback)
  {
    BOOST_REVERSE_FOREACH(const ScriptRunner::TracebackEntry& frame, *traceback)
    {
      if (frame.filename == wxT ("<string>"))
      {
	errorMarkerHandle = editor->MarkerAdd (frame.line-1, MARKER_Error);
	editor->Refresh ();
	return true;
      }
    }
    return false;
  }

  void ScriptEditorImpl::DisplayErrorLocation ()
  {
    if (errorMarkerHandle >= 0)
      editor->GotoLine (editor->MarkerLineFromHandle (errorMarkerHandle));
    GuiFrame::TabCallbackPtr callback (this->callback);
    if (callback) callback->Activate ();
    editor->SetFocus ();
  }
  
  void ScriptEditorImpl::SetupEditorStyle ()
  {
    editor->StyleClearAll();
    editor->SetLexer (wxSTC_LEX_PYTHON);
    for (int s = 0; s < wxSTC_STYLE_LASTPREDEFINED; s++)
    {
      editor->StyleSetFont (s, monospaceFont);
    }
    // Styles are taken from SciTE defaults
    editor->StyleSetSpec (wxSTC_P_DEFAULT, wxT ("fore:#808080"));
    editor->StyleSetSpec (wxSTC_P_COMMENTLINE, wxT ("fore:#007F00"));
    editor->StyleSetSpec (wxSTC_P_NUMBER, wxT ("fore:#007F7F"));
    editor->StyleSetSpec (wxSTC_P_STRING, wxT ("fore:#7F007F"));
    editor->StyleSetSpec (wxSTC_P_CHARACTER, wxT ("fore:#7F007F"));
    editor->StyleSetSpec (wxSTC_P_WORD, wxT ("bold,fore:#00007F"));
    editor->StyleSetSpec (wxSTC_P_TRIPLE, wxT ("fore:#7F0000"));
    editor->StyleSetSpec (wxSTC_P_TRIPLEDOUBLE, wxT ("fore:#7F0000"));
    editor->StyleSetSpec (wxSTC_P_CLASSNAME, wxT ("bold,fore:#0000FF"));
    editor->StyleSetSpec (wxSTC_P_DEFNAME, wxT ("bold,fore:#007F7F"));
    editor->StyleSetSpec (wxSTC_P_OPERATOR, wxT ("bold"));
    editor->StyleSetSpec (wxSTC_P_IDENTIFIER, wxT (""));
    editor->StyleSetSpec (wxSTC_P_COMMENTBLOCK, wxT ("fore:#7F7F7F"));
    editor->StyleSetSpec (wxSTC_P_STRINGEOL, wxT ("fore:#000000,back:#E0C0E0,eol"));
    editor->StyleSetSpec (wxSTC_P_WORD2, wxT ("fore:#203040,back:#FFFFC0"));
    editor->StyleSetSpec (wxSTC_P_DECORATOR, wxT ("fore:#805000"));
    editor->StyleSetSpec (wxSTC_STYLE_BRACELIGHT, wxT ("fore:#0000FF,bold"));
    editor->StyleSetSpec (wxSTC_STYLE_BRACEBAD, wxT ("fore:#FF0000,bold"));
    
    editor->SetKeyWords (0, wxT (
			 "and "
			 "as "
			 "assert "
			 "break "
			 "class "
			 "continue "
			 "def "
			 "del "
			 "elif "
			 "else "
			 "except "
			 "exec "
			 "finally "
			 "for "
			 "from "
			 "global "
			 "if "
			 "import "
			 "in "
			 "is "
			 "lambda "
			 "not "
			 "or "
			 "pass "
			 "print "
			 "raise "
			 "return "
			 "try "
			 "while "
			 "with "
			 "yield"
			 ));
			 
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDER,        wxSTC_MARK_ARROW,     wxT("BLACK"), wxT("BLACK"));
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDEROPEN,    wxSTC_MARK_ARROWDOWN, wxT("BLACK"), wxT("BLACK"));
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDERSUB,     wxSTC_MARK_VLINE,     wxT("#808080"), wxT("#808080"));
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDEREND,     wxSTC_MARK_ARROW,     wxT("BLACK"), wxT("BLACK"));
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDEROPENMID, wxSTC_MARK_ARROWDOWN, wxT("BLACK"), wxT("BLACK"));
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDERMIDTAIL, wxSTC_MARK_LCORNER, wxT("#808080"), wxT("#808080"));
    editor->MarkerDefine (wxSTC_MARKNUM_FOLDERTAIL,    wxSTC_MARK_LCORNER, wxT("#808080"), wxT("#808080"));
    
    editor->SetMarginType (marginLineNumber, wxSTC_MARGIN_NUMBER);
    
    editor->SetMarginType (marginFolding, wxSTC_MARGIN_SYMBOL);
    editor->SetMarginMask (marginFolding, wxSTC_MASK_FOLDERS);
    editor->SetMarginWidth (marginFolding, 16);
    editor->SetMarginSensitive (marginFolding, true);
    editor->SetProperty (wxT("fold"), wxT("1"));
    editor->SetProperty (wxT("fold.comment"), wxT("1"));
    editor->SetProperty (wxT("fold.compact"), wxT("1"));
    editor->SetProperty (wxT("fold.comment.python"), wxT("1"));
    editor->SetProperty (wxT("fold.quotes.python"), wxT("1"));
    editor->SetFoldFlags (wxSTC_FOLDFLAG_LINEBEFORE_CONTRACTED |
			  wxSTC_FOLDFLAG_LINEAFTER_CONTRACTED);
			  
    editor->SetIndent (2);
    
    editor->MarkerDefine (MARKER_Error, wxSTC_MARK_BACKGROUND, wxNullColour, wxColour (255, 192, 192));
  }
  
  void ScriptEditorImpl::UpdateLineNumberMargin()
  {
    int digitsNeeded = std::max (defaultLineNumberDigits,
				 int (ceil (log (editor->GetLineCount()) / log (10))));
    if (digitsNeeded == currentMarginDigits) return;
    
    editor->SetMarginWidth (marginLineNumber,
			    editor->TextWidth (wxSTC_STYLE_LINENUMBER, wxT ("0")) * (digitsNeeded + 0.5));
    
    currentMarginDigits = digitsNeeded;
  }

  wxString ScriptEditorImpl::GetDisplayFilename () const
  {
    if (fullPath.IsEmpty()) return wxT ("Unnamed Script");
    wxFileName fn (fullPath);
    return fn.GetFullName();
  }

  void ScriptEditorImpl::UpdateCaption ()
  {
    nutogui::GuiFrame::TabCallbackPtr callback (this->callback);
    if (!callback) return;
    
    wxString caption;
    if (isChanged) caption.Append (wxT ("*"));
    caption.Append (GetDisplayFilename ());
    callback->SetCaption (caption);
  }

  bool ScriptEditorImpl::DoFileOpen (const wxString& fullPath)
  {
    // @@@ Needed for WX 2.9
    editor->ClearAll();
    wxString textBackup (editor->GetText());
    
    if (!editor->LoadFile (fullPath))
    {
      editor->SetText (textBackup);
      uicommon::MessageDialog msgdlg;
      msgdlg.SetMessageHeading (wxString::Format (wxT ("Could not open %s."),
						  Quote::Single (fullPath).c_str()));
      msgdlg.SetIcon (wxICON_ERROR);
      msgdlg.OK (parentWindow);
      return false;
    }
    
    editor->SetSelection (0, 0); // @@@ Needed for WX 2.9
    
    wxFileName fn (fullPath);
    this->fullPath = fullPath;
    lastDir = fn.GetPath();
    fileHistory.AddFileToHistory (fullPath);
    SaveFilesConfig();
    isChanged = false;
    UpdateCaption ();
    return true;
  }
  
  bool ScriptEditorImpl::DoFileSave (bool saveas)
  {
    // Force 'Save as' when filename is empty
    saveas |= fullPath.IsEmpty();
    
    wxString saveFileName (fullPath);
    while (true)
    {
      if (saveas)
      {
	wxFileName saveFN (saveFileName);
	
	wxString saveDir (saveFN.GetPath());
	if (saveDir.IsEmpty()) saveDir = lastDir;
	wxFileDialog dlg (parentWindow, wxT("Save Script File"),
			  saveDir, 		// default dir
			  saveFN.GetFullName(), // default file
			  wxT ("Python Script|*.py"),
			  wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	
	if (dlg.ShowModal() != wxID_OK)
	  return false;
	saveFileName = dlg.GetPath ();
      }
      
      // Make a backup
      if (wxFileExists (saveFileName))
      {
	wxString backupName (saveFileName);
	backupName.Append (wxT ("~"));
	wxRenameFile (saveFileName, backupName, true);
      }
      bool fileOk = false;
      {
	wxTempFile tempFile (saveFileName);
	if (tempFile.IsOpened())
	{
	  // Always save files as UTF-8. @@@ Really?
	  wxCharBuffer textUTF8 (editor->GetText().ToUTF8());
	  size_t datalen = strlen (textUTF8.data());
	  fileOk = tempFile.Write (textUTF8.data(), datalen);
	  if (fileOk) fileOk = tempFile.Commit();
	}
      }
      if (fileOk)
      {
	wxFileName fn (saveFileName);
	fullPath = saveFileName;
	lastDir = fn.GetPath();
	fileHistory.AddFileToHistory (saveFileName);
	SaveFilesConfig();
	isChanged = false;
	UpdateCaption();
	return true;
      }
  
      // If saving failed, ask whether to save at some other location
      {
	uicommon::MessageDialog msgdlg;
	msgdlg.SetMessageHeading (wxString::Format (wxT ("Could not write to %s."),
						    Quote::Single (saveFileName).c_str()));
	msgdlg.SetMessageBody (wxT ("You may be able to save the file at another location."));
	msgdlg.SetQuestionDefault (wxT ("Would you like to choose another file name to save to?"));
	if (msgdlg.YesNo (parentWindow, wxID_YES, wxT ("Choose &another file name"), wxT ("&Don't save file")) == wxID_NO)
	  return false;
	
	saveas = true;
      }
      
      // Loop until we successfully save to a file or the user bailed out.
    }
  }
  
  bool ScriptEditorImpl::CanDiscardContents ()
  {
    if (!isChanged) return true;
    
    while (true)
    {
      uicommon::MessageDialog msgdlg;
      
      msgdlg.SetMessageHeading (wxString::Format (wxT ("%s was changed since it was last saved."),
						  Quote::Single (GetDisplayFilename()).c_str()));
      msgdlg.SetQuestionDefault (wxT ("Do you want to save the file now?"));
      int ret = msgdlg.YesNoCancel (parentWindow, wxID_CANCEL, wxID_SAVE, wxT ("&Don't save file"), wxID_CANCEL);
      switch (ret)
      {
      case wxID_CANCEL:	return false;
      case wxID_NO:	return true;
      }
      
      if (DoFileSave (false)) break;
    }
    
    return true;
  }
  
  void ScriptEditorImpl::SaveFilesConfig ()
  {
    wxConfigBase& config = *(wxConfigBase::Get ());
    config.SetPath (wxT ("ScriptEditor"));
    fileHistory.Save (config);
    config.Write (wxT ("lastDirectory"), lastDir);
    config.SetPath (wxT (".."));
  }
  
  void ScriptEditorImpl::OnFileNew (wxCommandEvent& event)
  {
    if (CanDiscardContents())
    {
      editor->ClearAll ();      
      fullPath.Clear();
      isChanged = false;
      UpdateCaption ();
    }
  }
  
  void ScriptEditorImpl::OnFileOpen (wxCommandEvent& event)
  {
    if (!CanDiscardContents()) return;
    
    wxFileDialog dlg (parentWindow, wxT("Open Script File"),
		      lastDir, 		// default dir
		      wxEmptyString, 	// default file
		      wxT ("Python Script|*.py"),
		      wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (dlg.ShowModal() == wxID_OK)
    {
      DoFileOpen (dlg.GetPath());
    }
  }

  void ScriptEditorImpl::OnFileHistory (wxCommandEvent& event)
  {
    if (!CanDiscardContents()) return;
    
    size_t idx = event.GetId()-wxID_FILE1;
    if (!DoFileOpen (fileHistory.GetHistoryFile (idx)))
    {
      fileHistory.RemoveFileFromHistory (idx);
      SaveFilesConfig();
    }
  }

  void ScriptEditorImpl::OnFileSave (wxCommandEvent& event)
  {
    DoFileSave (false);
  }
  
  void ScriptEditorImpl::OnFileSaveAs (wxCommandEvent& event)
  {
    DoFileSave (true);
  }
  
  void ScriptEditorImpl::OnArtificialDropDown (wxCommandEvent& event)
  {
    switch (event.GetId())
    {
    case ID_DropDownOpen:
      DoDropDown (wxID_OPEN);
      break;
    }
  }
    
  void ScriptEditorImpl::OnToolbarDropdown (wxAuiToolBarEvent& event)
  {
    if (!event.IsDropDownClicked())
    {
      return;
    }
    
    DoDropDown (event.GetToolId());
  }

  void ScriptEditorImpl::DoDropDown (int toolID)
  {
    wxPoint popupPos (toolbar->GetToolRect (toolID).GetBottomLeft());
    popupPos = toolbar->ClientToScreen (popupPos);
    popupPos = panel->ScreenToClient (popupPos);
    
    wxMenu* menu = nullptr;
    switch (toolID)
    {
    case wxID_OPEN:
      menu = &openPopup;
      break;
    case wxID_SAVE:
      menu = &savePopup;
      break;
    }
    
    if (menu)
      panel->PopupMenu (menu, popupPos);
  }
  
  void ScriptEditorImpl::OnCopy (wxCommandEvent& event)
  {
    editor->Copy ();
  }
  
  void ScriptEditorImpl::OnCut (wxCommandEvent& event)
  {
    editor->Cut ();
  }
  
  void ScriptEditorImpl::OnPaste (wxCommandEvent& event)
  {
    editor->Paste ();
  }
  
  void ScriptEditorImpl::OnUndo (wxCommandEvent& event)
  {
    editor->Undo ();
  }
  
  void ScriptEditorImpl::OnRedo (wxCommandEvent& event)
  {
    editor->Redo ();
  }
  
  void ScriptEditorImpl::OnUpdateUI (wxUpdateUIEvent& event)
  {
    switch (event.GetId())
    {
    case wxID_COPY:
    case wxID_CUT:
      {
	bool hasSelection = (editor->GetSelectionEnd() - editor->GetSelectionStart()) > 0;
	event.Enable (hasSelection);
      }
      break;
    case wxID_PASTE:
      event.Enable (editor->CanPaste());
      break;
    case wxID_UNDO:
      event.Enable (editor->CanUndo());
      break;
    case wxID_REDO:
      event.Enable (editor->CanRedo());
      break;
    }
  }

  void ScriptEditorImpl::OnEditorModified (wxStyledTextEvent& event)
  {
    bool newChangedState = editor->GetModify();
    if (isChanged != newChangedState)
    {
      isChanged = newChangedState;
      UpdateCaption ();
    }
  }
  
  static bool IsBrace (wxChar ch)
  {
    return (ch > 0) && (ch <= 0xff) && (strchr ("()[]{}<>", ch));
  }
  
  void ScriptEditorImpl::OnEditorUpdateUI (wxStyledTextEvent& event)
  {
    UpdateStatusBar();
    UpdateLineNumberMargin ();
    
    int currentPos = editor->GetCurrentPos();
    int curLinePos;
    wxCharBuffer curLine (editor->GetCurLineRaw (&curLinePos));
    bool brace = false;
    int braceOfs = 0;
    if (curLine.data())
    {
      // Check if current char (right of caret) is a brace
      int curChar = curLine[curLinePos];
      brace = IsBrace (curChar);
      // Current char is not a brace, but maybe the char left of the caret is
      if (!brace && (curLinePos > 0))
      {
	brace = IsBrace (curLine[curLinePos-1]);
	braceOfs = 1;
      }
    }
    if (brace)
    {
      int brace1 = currentPos - braceOfs;
      int brace2 = editor->BraceMatch (brace1);
      if (brace2 >= 0)
	editor->BraceHighlight (brace1, brace2);
      else
	editor->BraceBadLight (brace1);
    }
    else
      editor->BraceBadLight (-1);
    
    wxString highlightWord;
    int lexedStyle = editor->GetStyleAt (currentPos);
    if ((lexedStyle != wxSTC_P_IDENTIFIER) && (lexedStyle != wxSTC_P_WORD2))
    {
      // Maybe we're just off the end of a word... try a char to the left
      int col = editor->GetColumn (currentPos);
      if (col > 0)
      {
	currentPos--;
	lexedStyle = editor->GetStyleAt (currentPos);
      }
    }
    if ((lexedStyle == wxSTC_P_IDENTIFIER) || (lexedStyle == wxSTC_P_WORD2))
    {
      int wordStart = editor->WordStartPosition (currentPos, true);
      int wordEnd = editor->WordEndPosition (currentPos, true);
      highlightWord = editor->GetTextRange (wordStart, wordEnd);
    }
    newHighlightWord = highlightWord;
    /*if (highlightWord != currentHighlightWord)
    {
      currentHighlightWord = highlightWord;
      editor->SetKeyWords (1, highlightWord);
      editor->Colourise (0, editor->GetLength());
    }*/
  }
  
  void ScriptEditorImpl::UpdateStatusBar ()
  {
    positionIndicator->Freeze ();
    wxString statusText;
    if (editor->GetOvertype())
      statusText.Append (wxT ("OVR "));
    statusText.Append (wxString::Format (wxT ("Ln %d, Col %d"),
					 editor->GetCurrentLine() + 1,
					 editor->GetColumn (editor->GetCurrentPos()) + 1));
    positionIndicator->SetLabel (statusText);
    
    // Right align is broken on WX 2.8, so do it manually
    wxRect posIndicatorRect;
    statusBar->GetFieldRect (0, posIndicatorRect);
    posIndicatorRect.x = posIndicatorRect.GetRight() - positionIndicator->GetSize().x;
    positionIndicator->SetSize (posIndicatorRect);
    
    positionIndicator->Thaw ();
  }

  void ScriptEditorImpl::OnEditorMarginClick (wxStyledTextEvent& event)
  {
    if (event.GetMargin() == marginFolding)
    {
      int line = editor->LineFromPosition (event.GetPosition());
      int level = editor->GetFoldLevel (line);
      if ((level & wxSTC_FOLDLEVELHEADERFLAG) > 0)
	editor->ToggleFold (line);
    }
  }
  
  void ScriptEditorImpl::OnEditorCharAdded (wxStyledTextEvent& event)
  {
    int ch = event.GetKey();
    int currentLine = editor->GetCurrentLine();
    if (ch == '\n')
    {
	int lineIndent = 0;
	if (currentLine > 0)
	{
	  lineIndent = editor->GetLineIndentation (currentLine - 1);
	}
	if (lineIndent > 0)
	{
	  editor->SetUndoCollection (false);
	  editor->SetLineIndentation (currentLine, lineIndent);
	  editor->GotoPos (editor->PositionFromLine (currentLine) + lineIndent);
	  editor->SetUndoCollection (true);
	}
    }
  }
  
  void ScriptEditorImpl::OnSize (wxSizeEvent& event)
  {
    UpdateStatusBar ();
    
    event.Skip ();
  }
  
  void ScriptEditorImpl::OnIdle (wxIdleEvent& event)
  {
    if (newHighlightWord != currentHighlightWord)
    {
      currentHighlightWord = newHighlightWord;
      editor->SetKeyWords (1, newHighlightWord);
      editor->Colourise (0, editor->GetLength());
    }
    event.Skip();
  }
} // namespace nutogui
