%{
#include "mechanics/constitutive/damageLaws/DamageLaw.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/damageLaws/DamageLawLinear.h"
#include "mechanics/constitutive/damageLaws/DamageLawNoSoftening.h"
#include "mechanics/constitutive/damageLaws/DamageLawHermite.h"
%}

%include <std_shared_ptr.i>

%shared_ptr(NuTo::Constitutive::DamageLaw)
%shared_ptr(NuTo::Constitutive::DamageLawExponential)
%shared_ptr(NuTo::Constitutive::DamageLawLinear)
%shared_ptr(NuTo::Constitutive::DamageLawNoSoftening)
%shared_ptr(NuTo::Constitutive::DamageLawHermite)


%include "mechanics/constitutive/damageLaws/DamageLaw.h"
%include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
%include "mechanics/constitutive/damageLaws/DamageLawLinear.h"
%include "mechanics/constitutive/damageLaws/DamageLawNoSoftening.h"
%include "mechanics/constitutive/damageLaws/DamageLawHermite.h"
