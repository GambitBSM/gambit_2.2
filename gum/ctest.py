#!/usr/bin/env python

from lib.libfr import *

options = FROptions("feynrules", "SingletDM", "DiagonalCKM")

partlist = VectorOfParticles()
paramlist = VectorOfParameters()

all_feynrules(options, partlist, paramlist)

print partlist
print paramlist
