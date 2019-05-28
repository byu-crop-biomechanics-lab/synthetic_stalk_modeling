import numpy as np
from FEAclasses import *

Rind = Material("Rind",1.0e9,0.3)
Pith = Material("Pith",1.0e7,0.2)

Materials = [Rind,Pith]

# print(Rind.modulus)
