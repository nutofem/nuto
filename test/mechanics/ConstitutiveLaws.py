# -*- coding: utf-8 -*-
# $Id$

import sys
import nuto
import os
import math

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
createResult = False

#show the results on the screen
printResult = False

#system name and processor
system = sys.argv[1]+sys.argv[2]

#path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(sys.argv[3],"results",system,os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% start the real test file                                      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
structure=nuto.Structure(3)

Material1 = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
structure.ConstitutiveLawSetParameterDouble(Material1,"Youngs_Modulus", 20000)
structure.ConstitutiveLawSetParameterDouble(Material1,"Poissons_Ratio", 0.2)
structure.ConstitutiveLawSetParameterDouble(Material1,"Density", 0.5)
if(printResult):
    structure.ConstitutiveLawInfo(Material1, 0)
E = structure.ConstitutiveLawGetParameterDouble(Material1,"Youngs_Modulus")
if(E != 20000):
    print '[' + system,sys.argv[0] + '] : Young\'s modulus is not correct.'
    error = True
Nu = structure.ConstitutiveLawGetParameterDouble(Material1,"Poissons_Ratio")
if(Nu != 0.2):
    print '[' + system,sys.argv[0] + '] : Poisson\'s ratio is not correct.'
    error = True
Rho = structure.ConstitutiveLawGetParameterDouble(Material1,"Density")
if(Rho != 0.5):
    print '[' + system,sys.argv[0] + '] : Density is not correct.'
    error = True
        
structure.ConstitutiveLawDelete(Material1)
if(structure.GetNumConstitutiveLaws() != 0):
    print '[' + system,sys.argv[0] + '] : number of constitutive laws is not correct.'
    error = True
if(printResult):
    structure.ConstitutiveLawInfo(0)

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
