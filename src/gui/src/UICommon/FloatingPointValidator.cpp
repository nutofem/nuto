/**\file
 * wxValidator implementation for floating point input
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "uicommon/FloatingPointValidator.h"

#include <sstream>

namespace uicommon
{
  BEGIN_EVENT_TABLE(FloatingPointValidator, wxValidator)
    EVT_CHAR(FloatingPointValidator::OnChar)
  END_EVENT_TABLE()
  
  FloatingPointValidator::FloatingPointValidator (float& value, const std::locale& locale)
   : destination (value), locale (locale)
  {
  }

  bool FloatingPointValidator::ValidateOneChar (wxChar ch)
  {
    typedef std::numpunct<wxChar> NumPunctType;
    const NumPunctType& punctuation = std::use_facet<NumPunctType> (locale);
    return ((ch >= '0') && (ch <= '9'))
      || (ch == '+') || (ch == '-')
      || (ch == 'E') || (ch == 'e')
      || (ch == punctuation.decimal_point())
      || (ch == punctuation.thousands_sep());
  }
  
  class num_get_from_wxString : public std::num_get<wxChar, wxString::const_iterator>
  {
  public:
    num_get_from_wxString() : std::num_get<wxChar, wxString::const_iterator> (1) {}
  };

  bool FloatingPointValidator::TryConvert (const wxString& str, float& v)
  {
    std::ios_base::iostate err (std::ios_base::goodbit);
    std::istringstream formatting;
    formatting.imbue (locale);
    num_get_from_wxString getter;
    getter.get (str.begin(), str.end(), formatting, err, v);
    return ((err & std::ios_base::failbit) == 0);
  }

  wxObject* FloatingPointValidator::Clone () const
  {
    return new FloatingPointValidator (destination, locale);
  }
  
  bool FloatingPointValidator::TransferFromWindow ()
  {
    wxTextCtrl* textCtrl = wxDynamicCast (GetWindow (), wxTextCtrl);
    wxCHECK2_MSG(textCtrl, return false, "FloatingPointValidator expects wxTextCtrl windows");

    float v;
    if (!TryConvert (textCtrl->GetValue(), v)) return false;
    destination = v;
    return true;
  }
  
  bool FloatingPointValidator::TransferToWindow ()
  {
    wxTextCtrl* textCtrl = wxDynamicCast (GetWindow (), wxTextCtrl);
    wxCHECK2_MSG(textCtrl, return false, "FloatingPointValidator expects wxTextCtrl windows");
    
    std::stringstream outstream;
    outstream.imbue (locale);
    typedef std::num_put<char> NumPutter;
    std::use_facet<NumPutter> (locale).put (outstream, outstream, ' ', destination);
    if (outstream.fail()) return false;
    
    textCtrl->ChangeValue (wxString (outstream.str().c_str(), wxConvLocal));
    return true;
  }
  
  bool FloatingPointValidator::Validate (wxWindow* parent)
  {
    wxTextCtrl* textCtrl = wxDynamicCast (GetWindow (), wxTextCtrl);
    wxCHECK2_MSG(textCtrl, return false, "FloatingPointValidator expects wxTextCtrl windows");

    float v;
    return TryConvert (textCtrl->GetValue(), v);
  }

  void FloatingPointValidator::OnChar (wxKeyEvent& event)
  {
    int keyCode = event.GetKeyCode();
    int keyChar;
  #if defined(wxHAS_UNICODE) && wxHAS_UNICODE
    keyChar = event.GetUnicodeKey();
  #else
    keyChar = keyCode;
  #endif

    if ((keyChar >= WXK_SPACE) && (keyChar < WXK_START) && (keyChar != WXK_DELETE))
    {
      if (!ValidateOneChar (keyChar))
      {
        if (!wxValidator::IsSilent())
	  wxBell();
	return;
      }
    }
    
    event.Skip();
  }

} // namespace uicommon
