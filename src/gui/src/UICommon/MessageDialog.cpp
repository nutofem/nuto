#include "common.h"

#include "uicommon/MessageDialog.h"

#include <wx/wx.h>

namespace uicommon
{
  void MessageDialog::SetupStrings (wxString& caption,
				    wxString& initMessageStr)
  {
    caption = this->caption;
  #if wxABI_VERSION >= 20900
    if (caption.IsEmpty ()) caption = wxTheApp->GetAppDisplayName ();
  #else
    if (caption.IsEmpty ()) caption = wxTheApp->GetAppName ();
  #endif
    
    initMessageStr = heading;
  #if wxABI_VERSION < 20900
    if (!body.IsEmpty ())
    {
      initMessageStr.Append (wxT ("\n"));
      initMessageStr.Append (body);
    }
    if (!questionDefault.IsEmpty())
    {
      initMessageStr.Append (wxT ("\n"));
      initMessageStr.Append (questionDefault);
    }
  #endif
  }
  
  void MessageDialog::SetupExtMessage (wxMessageDialog& msgdlg,
				       bool labelsSet)
  {
  #if wxABI_VERSION >= 20900
    wxString extmsg (body);
    wxString questionStr;
    if (labelsSet)
      questionStr = !questionCustom.IsEmpty() ? questionCustom : questionDefault;
    else
      questionStr = questionDefault;
    if (!questionStr.IsEmpty())
    {
      if (!extmsg.IsEmpty ()) extmsg.Append (wxT ("\n"));
      extmsg.Append (questionStr);
    }
    msgdlg.SetExtendedMessage (extmsg);
  #endif
  }
  
  #if wxABI_VERSION >= 20900
    #define SET_LABELS(Expr)	Expr
  #else
    #define SET_LABELS(Expr)	false
  #endif
  
  void MessageDialog::OK (wxWindow* parent,
			  const ButtonLabel& customLabelOK)
  {
    wxString caption, initMessageStr;
    SetupStrings (caption, initMessageStr);
    
    long style = iconID | wxOK;
    wxMessageDialog msgdlg (parent, initMessageStr, caption, style);
    SetupExtMessage (msgdlg,
		     SET_LABELS (msgdlg.SetOKLabel (customLabelOK)));
    msgdlg.ShowModal ();
  }
  
  int MessageDialog::YesNo (wxWindow* parent,
			    int defaultLabel,
			    const ButtonLabel& customLabelYes,
			    const ButtonLabel& customLabelNo)
  {
    wxString caption, initMessageStr;
    SetupStrings (caption, initMessageStr);
    
    long style = iconID | wxYES_NO;
    switch (defaultLabel)
    {
    case wxID_YES:	style |= wxYES_DEFAULT;
    case wxID_NO:	style |= wxNO_DEFAULT;
    }
    wxMessageDialog msgdlg (parent, initMessageStr, caption, style);
    SetupExtMessage (msgdlg,
		     SET_LABELS (msgdlg.SetYesNoLabels (customLabelYes, customLabelNo)));
    return msgdlg.ShowModal ();
  }
  
  int MessageDialog::YesNoCancel (wxWindow* parent,
				  int defaultLabel,
				  const ButtonLabel& customLabelYes,
				  const ButtonLabel& customLabelNo,
				  const ButtonLabel& customLabelCancel)
  {
    wxString caption, initMessageStr;
    SetupStrings (caption, initMessageStr);
    
    long style = iconID | wxYES_NO | wxCANCEL;
    switch (defaultLabel)
    {
    case wxID_YES:	style |= wxYES_DEFAULT;
    case wxID_NO:	style |= wxNO_DEFAULT;
  #if wxABI_VERSION >= 20900
    case wxID_CANCEL:	style |= wxCANCEL_DEFAULT;
  #endif
    }
    wxMessageDialog msgdlg (parent, initMessageStr, caption, style);
    SetupExtMessage (msgdlg,
		     SET_LABELS (msgdlg.SetYesNoCancelLabels (customLabelYes, customLabelNo, customLabelCancel)));
    return msgdlg.ShowModal ();
  }
} // namespace uicommon
