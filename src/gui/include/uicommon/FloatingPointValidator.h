/**\file
 * wxValidator implementation for floating point input
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_FLOATINGPOINTVALIDATOR_H__
#define __UICOMMON_FLOATINGPOINTVALIDATOR_H__

#include "export.h"

#include <wx/wx.h>

#include <locale>

namespace uicommon
{
  /**
   * wxValidator implementation for floating point input in text controls.
   * Like standard WX validators, it reads the value of a text control and converts
   * it to a float value, or vice versa.
   */
  class UICOMMON_EXPORTED FloatingPointValidator : public wxValidator
  {
    float& destination;
    const std::locale& locale;
    
    bool ValidateOneChar (wxChar ch);
    bool TryConvert (const wxString& str, float& v);
    void OnChar (wxKeyEvent& event);
  public:
    /**
     * Construct.
     * \param value Reference to float variable that is read from when transferring to
     *   the control resp. written to when transferring from the control.
     * \param locale C++ locale to use when converting floats to or from string.
     */
    FloatingPointValidator (float& value, const std::locale& locale);
    
    wxObject* Clone () const;
    bool TransferFromWindow ();
    bool TransferToWindow ();
    bool Validate (wxWindow* parent);

    DECLARE_EVENT_TABLE()
  };
} // namespace uicommon

#endif // __UICOMMON_FLOATINGPOINTVALIDATOR_H__
