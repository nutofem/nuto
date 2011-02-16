#include <vtkRenderWindowInteractor.h>

class Interactor : public vtkRenderWindowInteractor
{
public:
  using vtkRenderWindowInteractor::InternalCreateTimer;
};

int main()
{
  Interactor* rwi = 0;
  rwi->InternalCreateTimer (0, 0, 0);
  return 1;
}
