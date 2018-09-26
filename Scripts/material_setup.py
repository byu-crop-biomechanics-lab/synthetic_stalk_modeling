# Save by ryanal on 2018_09_26-16.27.20; build 2018 2017_11_07-10.21.41 127140
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

class Material:
    def __init__(self,name,modulus,poisson):
        self.name = name
        self.modulus = modulus
        self.poisson = poisson

Rind = Material("Rind",1.0e9,0.3)
Pith = Material("Pith",1.0e7,0.2256)

mp = [Rind,Pith]

for material in mp:
    mdb.models['Model-1'].Material(name=material.name)
    mdb.models['Model-1'].materials[material.name].Elastic(table=((material.modulus,
                                                        material.poisson), ))

# Save by ryanal on 2018_09_26-16.28.29; build 2018 2017_11_07-10.21.41 127140
