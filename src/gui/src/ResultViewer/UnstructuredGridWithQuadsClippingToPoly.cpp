/**\file
 * vtkUnstructuredGrid subclass that contains QuadClippingToPoly instances for quads.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "UnstructuredGridWithQuadsClippingToPoly.h"

#include "QuadClippingToPoly.h"

#include <vtkGenericCell.h>
#include <vtkUnsignedCharArray.h>

namespace nutogui
{
  /*
    'Customized' unstructured grid clipping.
    The normal vtkQuad can clip to multiple tris or quads. However, in edge
    drawing mode, this results in unwanted additional edges.
    To this end, QuadClippingToPoly is used, which clips to one poly - thus,
    no additional edges.
    The 'Hack*' classes below are just means to have QuadClippingToPoly instead
    of vtkQuad used in the clipping.
    
    (Alternatively, a custom vtkQuad object factory could've been used. However,
    swapping out the vtkQuad implementation globally may have other side effects...
    so use the roundabout, but 'contained' way below.)
   */
  
  class HackGenericCell : public vtkGenericCell
  {
  public:
    /* Overridden to output QuadClippingToPoly instead of vtkQuad.
       *Not* a virtual method!... */
    void SetCellType (int cellType);
  };
  
  vtkCxxRevisionMacro(UnstructuredGridWithQuadsClippingToPoly, "1")
  
  void UnstructuredGridWithQuadsClippingToPoly::GetCell (vtkIdType cellId, vtkGenericCell *cell)
  {
    /* vtkGenericCell::SetCellType() is *not* virtual, so manually call
       our 'customized' variant */
    static_cast<HackGenericCell*> (cell)->SetCellType (Types->GetValue(cellId));
    
    /* Superclass GetCell() will also set the cell type,
       but that will not change the 'Cell' member. */
    vtkUnstructuredGrid::GetCell (cellId, cell);
  }
  
  void HackGenericCell::SetCellType (int cellType)
  {
    if ( this->Cell->GetCellType() != cellType )
    {
      this->Points->UnRegister(this);
      this->PointIds->UnRegister(this);
      this->PointIds = NULL;
      this->Cell->Delete();

      vtkCell* cell;
      if (cellType == VTK_QUAD)
	// Our own quad
	cell = QuadClippingToPoly::New();
      else
	// Default instantiation.
	cell = vtkGenericCell::InstantiateCell(cellType);

      this->Cell = cell;
      this->Points = this->Cell->Points;
      this->Points->Register(this);
      this->PointIds = this->Cell->PointIds;
      this->PointIds->Register(this);
    }//need to change cell type
  }

} // namespace nutogui
