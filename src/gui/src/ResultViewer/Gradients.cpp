/**\file
 * List of gradients for result viewer view panel.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "Gradients.h"

#include "uicommon/BitmapGraphicsContext.h"

#include <wx/graphics.h>

#include <vtkColorTransferFunction.h>
#include <vtkLookupTable.h>

namespace nutogui
{
  static vtkSmartPointer<vtkLookupTable> ColorTransferFunctionToLookupTable (vtkColorTransferFunction* ctf)
  {
    size_t nColors = 256;
    vtkSmartPointer<vtkLookupTable> colors (vtkSmartPointer<vtkLookupTable>::New());
    colors->SetNumberOfTableValues (nColors);
    for (size_t i = 0; i < nColors; i++)
    {
      double rgb[3];
      ctf->GetColor (float (i)/float (nColors-1), rgb);
      colors->SetTableValue (i, rgb[0], rgb[1], rgb[2]);
    }
    return colors;
  }
  
  ResultViewerImpl::View::Gradients::Gradients ()
  {
    // Set up default gradients
    vtkSmartPointer<vtkLookupTable> colors;
    vtkSmartPointer<vtkColorTransferFunction> ctf;
    {
      colors = vtkSmartPointer<vtkLookupTable>::New(); 
      colors->SetHueRange (0.667, 0.0);
      colors->Build ();
      
      gradients.push_back (std::make_pair (wxT("Rainbow"), colors));
    }
    {
      ctf = vtkSmartPointer<vtkColorTransferFunction>::New(); 
      ctf->AddRGBPoint (0, 0, 0, 0);
      ctf->AddRGBPoint (1, 1, 1, 1);
      colors = ColorTransferFunctionToLookupTable (ctf);
      
      gradients.push_back (std::make_pair (wxT("Grayscale"), colors));
    }
    // "pm3d" gradient as e.g. found in gnuplot
    {
      size_t nColors = 256;
      vtkSmartPointer<vtkLookupTable> colors (vtkSmartPointer<vtkLookupTable>::New());
      colors->SetNumberOfTableValues (nColors);
      for (size_t i = 0; i < nColors; i++)
      {
	double x = i / double (nColors-1);
	colors->SetTableValue (i, sqrt (x), x*x*x, sin (x*M_PI));
      }
      
      gradients.push_back (std::make_pair (wxT("pm3d"), colors));
    }
  }
  
  vtkLookupTable* ResultViewerImpl::View::Gradients::GetGradientColors (size_t n) const
  {
    return gradients[n].second;
  }
  
  wxBitmap ResultViewerImpl::View::Gradients::RenderGradient (size_t n, const wxSize& bitmapSize) const
  {
    vtkScalarsToColors* colors = gradients[n].second;
    
    // Direction of gradient *strip*: vertical gradient -> horizontal strips, horizontal gradient -> vertical strips!
    bool stripsHorizontal = true;
    wxBitmap bmp (bitmapSize.x, bitmapSize.y, 32);
    {
      uicommon::BitmapGraphicsContext gc (bmp);

      gc->SetPen (*wxTRANSPARENT_PEN);
      int numStrips = (stripsHorizontal ? bitmapSize.y : bitmapSize.x) - 2;
      int stripW = stripsHorizontal ? bitmapSize.x - 2 : 1;
      int stripH = stripsHorizontal ? 1 : bitmapSize.y - 2;
      for (int strip = 0; strip < numStrips; strip++)
      {
	int x, y;
	if (stripsHorizontal)
	{
	  x = 1;
	  y = strip + 1;
	}
	else
	{
	  x = strip + 1;
	  y = 1;
	}
	
	double rgb[3];
	double v = strip / double (numStrips-1);
	if (stripsHorizontal) v = 1-v;
	colors->GetColor (v, rgb);
	gc->SetBrush (wxColour (rgb[0]*255, rgb[1]*255, rgb[2]*255));
	gc->DrawRectangle (x, y, stripW, stripH);
      }
      gc->SetBrush (*wxTRANSPARENT_BRUSH);
      {
	gc->SetPen (wxColour (128, 128, 128));
	gc->DrawRoundedRectangle (0, 0, bitmapSize.x-1, bitmapSize.y-1, 1); // ??? -1 -- WX off by one?
      }
    }
    return bmp;
  }
} // namespace nutogui
