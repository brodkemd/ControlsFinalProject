from deflectioncalc import deflection_calc
import numpy as np

df = deflection_calc()

dels = df.calc_def(np.array([0,0,0]),np.array([0,0,0,5.0,0,150.0]))

