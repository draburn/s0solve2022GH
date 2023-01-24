import numpy as np

matD = np.zeros(( 4, 3 ))
matX = np.array([[ 11, 12, 13 ],
                 [ 21, 22, 23 ],
                 [ 31, 32, 33 ],
                 [ 41, 42, 43 ]])
vecXAnchor = np.array([ 10, 20, 30, 40 ])
#vecXAnchor = np.array([[10], [20], [30], [40] ])
n = 3
matD[:,0:n] = matX[:,0:n] - np.reshape( vecXAnchor, (4,1) )
print( matD )


from scipy import linalg
matV = linalg.orth( matX )
print( matX )
print( matV )
print( matV.T @ matX )

# Consider using overwrite_a = True with qr()
matQ, matR = linalg.qr( matX, mode='economic' )
print( matX )
print( matQ )
print( matR )
print( matQ.T @ matX )
