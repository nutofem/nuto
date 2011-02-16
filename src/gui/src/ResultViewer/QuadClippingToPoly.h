/**\file
 * Subclass of vtkQuad with different clipping output.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_QUADCLIPPINGTOPOLY_H__
#define __NUTOGUI_QUADCLIPPINGTOPOLY_H__

#include <vtkQuad.h>

namespace nutogui
{
  /**
   * A vtkQuad which outputs one polygon when clipped.
   * The regulat vtkQuad outputs up to three triangles or quads.
   * This variant outputs a single polygon (but with up to 6 vertices).
   * Useful when the 'implicit' triangulation done by vtkQuad is undesireable.
   */
  class QuadClippingToPoly : public vtkQuad
  {
  public:
    static QuadClippingToPoly* New ()
    { return new QuadClippingToPoly; }
    
    void Clip(double value, vtkDataArray *cellScalars, 
	      vtkIncrementalPointLocator *locator, vtkCellArray *polys,
	      vtkPointData *inPd, vtkPointData *outPd,
	      vtkCellData *inCd, vtkIdType cellId, vtkCellData *outCd,
	      int insideOut);
  };
} // namespace nutogui

#endif // __NUTOGUI_QUADCLIPPINGTOPOLY_H__
