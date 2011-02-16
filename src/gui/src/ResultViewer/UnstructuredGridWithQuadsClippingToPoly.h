/**\file
 * vtkUnstructuredGrid subclass that contains QuadClippingToPoly instances for quads.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_UNSTRUCTUREDGRIDWITHQUADSCLIPPINGTOPOLY_H__
#define __NUTOGUI_UNSTRUCTUREDGRIDWITHQUADSCLIPPINGTOPOLY_H__

#include "vtkUnstructuredGrid.h"

namespace nutogui
{
  /**
   * vtkUnstructuredGrid subclass that contains QuadClippingToPoly instances for quads.
   * Useful for clipping.
   */
  class UnstructuredGridWithQuadsClippingToPoly : public vtkUnstructuredGrid
  {
  public:
    static UnstructuredGridWithQuadsClippingToPoly* New()
    { return new UnstructuredGridWithQuadsClippingToPoly; }

    vtkTypeRevisionMacro(UnstructuredGridWithQuadsClippingToPoly,vtkUnstructuredGrid);
    
    // Overridden to output QuadClippingToPoly instead of vtkQuad.
    void GetCell (vtkIdType cellId, vtkGenericCell *cell);
  };
} // namespace nutogui

#endif // __NUTOGUI_UNSTRUCTUREDGRIDWITHQUADSCLIPPINGTOPOLY_H__
  