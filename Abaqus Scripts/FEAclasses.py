import numpy as np

class Material:
    def __init__(self,name,modulus,poisson):
        self.name = name
        self.modulus = modulus
        self.poisson = poisson

