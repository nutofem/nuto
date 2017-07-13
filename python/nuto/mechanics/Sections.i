// ignore << operator, because SWIG issues a warning otherwise and we're wrapping it manually (see %extend below)
%ignore operator<<;

%{
#include "mechanics/sections/SectionFibreMatrixBond.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/sections/SectionVariableTruss.h"
%}

%include <std_shared_ptr.i>

%shared_ptr(NuTo::Section)
%shared_ptr(NuTo::SectionTruss)
%shared_ptr(NuTo::SectionPlane)
%shared_ptr(NuTo::SectionFibreMatrixBond)
%shared_ptr(NuTo::SectionVariableTruss)

%include "mechanics/sections/Section.h"
%include "mechanics/sections/SectionFibreMatrixBond.h"
%include "mechanics/sections/SectionPlane.h"
%include "mechanics/sections/SectionTruss.h"
%include "mechanics/sections/SectionVariableTruss.h"
