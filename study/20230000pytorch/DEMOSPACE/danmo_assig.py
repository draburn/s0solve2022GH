import numpy as np

matJ = np.zeros( (2,2) )
vecD = np.zeros( 2 )
vecG0 = np.zeros( 2 )
vecGA = np.zeros( 2 )
vecGB = np.zeros( 2 )

vecG1 = vecG0
vecG2 = matJ @ vecD
vecG3 = 1*vecG0
vecGA = vecG0
vecGB[:] = vecG0[:]

vecG1[0] = 1
vecG2[0] = 2
vecG3[0] = 3

matJ
vecD
vecG0
vecG1
vecG2
vecG3
vecGA
vecGB
