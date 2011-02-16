/**\file
 * Tab to display output from a script
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "ScriptOutputTab.h"

#include <wx/richtext/richtextctrl.h>
#include <boost/weak_ptr.hpp>

ScriptOutputTab::ScriptOutputTab () : lastStyle (0)
{
}

wxWindow* ScriptOutputTab::CreateContents (const nutogui::GuiFrame::TabCallbackWeakPtr& callback_,
					   wxWindow* parentWindow)
{
  nutogui::GuiFrame::TabCallbackPtr callback (callback_);
  callback->SetCaption (wxT ("Script Output"));
  
  wxPanel* panel = new wxPanel (parentWindow);
  wxSizer* sizer = new wxBoxSizer (wxVERTICAL);

  textCtrl = new wxTextCtrl (panel, wxID_ANY, wxEmptyString,
			     wxDefaultPosition, wxDefaultSize,
			     wxTE_READONLY | wxTE_MULTILINE | wxTE_RICH2 | wxBORDER_NONE);
  sizer->Add (textCtrl, 1, wxEXPAND);
  
  panel->SetSizer (sizer);
  
  {
    wxFont defFont (wxSystemSettings::GetFont (wxSYS_DEFAULT_GUI_FONT));
    wxFont monospaceFont (wxFont (defFont.GetPointSize(), wxFONTFAMILY_MODERN, defFont.GetStyle(), defFont.GetWeight()));
    textCtrl->SetFont (monospaceFont);
    textCtrl->SetBackgroundColour (wxSystemSettings::GetColour (wxSYS_COLOUR_BTNFACE));
    wxTextAttr textAttr;
    textAttr.SetFlags (wxTEXT_ATTR_TEXT_COLOUR);
    textAttr.SetTextColour (wxSystemSettings::GetColour (wxSYS_COLOUR_WINDOWTEXT));
    textCtrl->SetDefaultStyle (textAttr);
  }
  
  return panel;
}

void ScriptOutputTab::AppendOutput (const wxString& output)
{
  AppendOutput (0, output);
}

void ScriptOutputTab::AppendErrorOutput (const wxString& output)
{
  AppendOutput (lineError, output);
}

void ScriptOutputTab::AppendNote (const wxString& output)
{
  AppendOutput (lineNote, output);
}

void ScriptOutputTab::AppendOutput (unsigned int style, const wxString& output)
{
  wxTextAttr textAttr;
  int attrFlags = 0;
  if ((style & lineError) != (lastStyle & lineError))
  {
    attrFlags |= wxTEXT_ATTR_TEXT_COLOUR;
    if (style & lineError)
      textAttr.SetTextColour (wxColour (128, 0, 0));
    else
      textAttr.SetTextColour (wxSystemSettings::GetColour (wxSYS_COLOUR_WINDOWTEXT));
  }
  if ((style & fontStyles) != (lastStyle & fontStyles))
  {
    wxFont font (textAttr.HasFont() ? textAttr.GetFont() : textCtrl->GetFont ());
    int fontFlags = 0;
    if ((style & lineNote) != (lastStyle & lineNote))
    {
      fontFlags |= wxTEXT_ATTR_FONT_ITALIC;
      if (style & lineNote)
	font.SetStyle (wxFONTSTYLE_ITALIC);
      else
	font.SetStyle (wxFONTSTYLE_NORMAL);
    }
    textAttr.SetFont (font, 0);
    attrFlags |= fontFlags;
  }
  if (attrFlags != 0)
  {
    textAttr.SetFlags (attrFlags);
    textCtrl->SetDefaultStyle (textAttr);
  }
  lastStyle = style;
    
  long oldSelFrom, oldSelTo;
  textCtrl->GetSelection (&oldSelFrom, &oldSelTo);
  long oldInsPt = textCtrl->GetInsertionPoint();
  long lastPos = textCtrl->GetLastPosition();
  textCtrl->AppendText (output);
  // If the insertion point was not at the end, restore it (so to not annoy the user)
  if (oldInsPt != lastPos)
    textCtrl->SetInsertionPoint (oldInsPt);
  // Restore selection
  if (oldSelTo > oldSelFrom)
    textCtrl->SetSelection (oldSelFrom, oldSelTo);
}
