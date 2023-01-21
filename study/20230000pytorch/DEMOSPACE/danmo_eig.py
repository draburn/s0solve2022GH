# https://numpy.org/doc/stable/reference/generated/numpy.linalg.eig.html

import numpy as np
from numpy import linalg as LA
w, v = LA.eig(np.diag((1,2,3)))
print( "w = ", w )
print( "v = ", v )
