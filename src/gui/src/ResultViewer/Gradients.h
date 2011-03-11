/**\file
 * List of gradients for result viewer view panel.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __GRADIENTS_H__
#define __GRADIENTS_H__

#include "View.h"

#include <wx/string.h>

#include <vtkSmartPointer.h>

class vtkLookupTable;

namespace nutogui
{
  /**
   * Manage list of gradients.
   */
  class ResultViewerImpl::View::Gradients
  {
    std::vector<std::pair<wxString, vtkSmartPointer<vtkLookupTable> > > gradients;
  public:
    Gradients ();
    
    /// Number of available gradients
    size_t GetGradientsNum () const
    { return gradients.size(); }
    /// Get name of a gradient
    const wxString& GetGradientName (size_t n) const
    { return gradients[n].first; }
    /// Get VTK colors for gradient.
    vtkLookupTable* GetGradientColors (size_t n) const;
    
    /// Render a gradient to a bitmap.
    wxBitmap RenderGradient (size_t n, const wxSize& size) const;
  };
} // namespace nutogui

#endif // __GRADIENTS_H__
