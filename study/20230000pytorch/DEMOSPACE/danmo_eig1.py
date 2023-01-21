import numpy as np
from numpy import linalg as LA
m = np.array( [[1, 1], [1,-1]] )
w, v = LA.eig(m)
print( "m = ", m )
print( "w = ", w )
print( "v = ", v )
