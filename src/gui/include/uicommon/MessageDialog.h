#ifndef __UICOMMON_MESSAGEDIALOG_H__
#define __UICOMMON_MESSAGEDIALOG_H__

#include "export.h"

#include <wx/wx.h>
#include <wx/msgdlg.h>

namespace uicommon
{
  /**
   * A wrapper for wxMessageDialog to allow a unified interface,
   * independent of whether WX 2.8 or 2.9 is used (the latter
   * offering much more features).
   */
  class UICOMMON_EXPORTED MessageDialog
  {
    long iconID;
    wxString caption;
    wxString heading;
    wxString body;
    
    wxString questionDefault;
    wxString questionCustom;
    
    void SetupStrings (wxString& caption, wxString& initMessageStr);
    void SetupExtMessage (wxMessageDialog& msgdlg, bool labelsSet);
  public:
  #if wxABI_VERSION >= 20900
    typedef wxMessageDialog::ButtonLabel ButtonLabel;
  #else
    struct ButtonLabel
    {
      ButtonLabel (int) {}
      ButtonLabel (const wxChar*) {}
    };
  #endif
  
    MessageDialog () : iconID (0) {}
    
    /// Set icon to use
    void SetIcon (long iconID)
    { this->iconID = iconID; }
    /// Set dialog caption (defaults to app name)
    void SetCaption (const wxString& caption)
    { this->caption = caption; }
    /// Set message heading string  (displayed first with some emphasis; should be short)
    void SetMessageHeading (const wxString& heading)
    { this->heading = heading; }
    /// Set message body string (displayed after heading; can be longer/explanatory)
    void SetMessageBody (const wxString& body)
    { this->body = body; }
    
    /// Set question string, used when only default button labels are available
    void SetQuestionDefault (const wxString& question)
    { this->questionDefault = question; }
    /// Set question string, used when custom button labels are available
    void SetQuestionCustom (const wxString& question)
    { this->questionCustom = question; }
    
    /// Display a notification box
    void OK (wxWindow* parent,
	     const ButtonLabel& customLabelOK = wxID_OK);
    /// Ask a "Yes"/"No" style question, with optional custom buttons
    int YesNo (wxWindow* parent,
	       int defaultLabel = wxID_YES,
	       const ButtonLabel& customLabelYes = wxID_YES,
	       const ButtonLabel& customLabelNo = wxID_NO);
    /// Ask a "Yes"/"No"/"Cancel" style question, with optional custom buttons
    int YesNoCancel (wxWindow* parent,
		     int defaultLabel = wxID_YES,
		     const ButtonLabel& customLabelYes = wxID_YES,
		     const ButtonLabel& customLabelNo = wxID_NO,
		     const ButtonLabel& customLabelCancel = wxID_CANCEL);
  };
} // namespace uicommon

#endif // __UICOMMON_MESSAGEDIALOG_H__
