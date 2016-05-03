import numpy as np
import scipy.stats as sp
x=np.rand(1000,1)
y=np.rand(1000,1)

print(sp.pearsonr(x, y))
