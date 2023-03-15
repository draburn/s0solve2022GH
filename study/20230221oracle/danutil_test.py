import danutil
from danutil import msg
import numpy as np
from numpy.linalg import norm

a = np.random.randn()
b = np.random.randn()
c = np.random.randn()
msg(f'a = {a}, b = {b}, c = {c}')
xlroq = danutil.linishrootOfQuad(a, b, c)
xLo = np.floor(-abs(xlroq))
xHi = np.ceil(abs(xlroq))
numVals = 1001
xVals = np.linspace(xLo, xHi, 1001)
fVals = (a*(xVals**2)) + (b*xVals) + c
flroq = a*(xlroq**2) + b*xlroq + c
msg(f'xlroq = {xlroq}, flroq = {flroq}')
gVals = b*xVals + c

import matplotlib.pyplot as plt
plt.plot( [xLo, xHi], [0.0, 0.0], '-', color='black' )
plt.plot( [0.0, 0.0], [np.floor(min(fVals)), np.ceil(max(fVals))], '-', color='black' )
plt.plot( 0.0, c, 's', markersize=10 )
plt.plot( xVals, gVals, '-' )
plt.plot( xVals, fVals, '.-' )
plt.plot( xlroq, flroq, '*', markersize=15 )
plt.grid(True)
plt.show()
