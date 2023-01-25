import numpy as np

matD = np.zeros(( 4, 3 ))
matX = np.array([[ 11, 12, 13 ],
                 [ 21, 22, 23 ],
                 [ 31, 32, 33 ],
                 [ 41, 42, 43 ]])
#matX = np.array([[ 1.1, 1.1, 1.1 ],
#                 [ 2.1, 2.1, 2.1 ],
#                 [ 3.1, 3.1+1E-8, 3.1 ],
#                 [ 4.1, 4.1, 4.2 ]])
vecXAnchor = np.array([ 10, 20, 30, 40 ])
#vecXAnchor = np.array([[10], [20], [30], [40] ])
n = 3
matD[:,0:n] = matX[:,0:n] - np.reshape( vecXAnchor, (4,1) )
print( matD )

print()
print( 'ortho()...' )
from scipy import linalg
#myNorm = linalg.norm(matX[:,-1])
#print( 'myNorm = ', myNorm )
#matV = linalg.orth( matX, rcond=None ) # Doesn't work!
matV = linalg.orth( matX )
print( 'matX =\n', matX )
print( 'matV =\n', matV )
print( 'V^T * X =\n', matV.T @ matX )

print()
print( 'Sans pivoting...' )
# Consider using overwrite_a = True with qr()
matQ, matR = linalg.qr( matX, mode='economic' )
print( 'matX =\n', matX )
print( 'matQ =\n', matQ )
print( 'matR =\n', matR )
print( 'Q^T * X =\n', matQ.T @ matX )

print()
print( 'With pivoting...' )
matQ, matR, matP = linalg.qr( matX, mode='economic', pivoting=True )
print( 'matX =\n', matX )
print( 'matQ =\n', matQ )
print( 'matR =\n', matR )
print( 'matP =\n', matP )
print( 'Q^T * X =\n', matQ.T @ matX )
print( 'Q^T * X[P] =\n', matQ.T @ ( (matX.T[matP]).T ) )
print( 'X[P] =\n', (matX.T[matP]).T )

print()
print( 'Two-pass QR...' )
# Consider using overwrite_a = True with qr()
matQ1, matR1 = linalg.qr( matX, mode='economic' )
print( 'matX =\n', matX )
print( 'matQ1 =\n', matQ1 )
print( 'matR1 =\n', matR1 )
#rvcKeep = np.abs(np.diag(matR1)) > 0.1
dropThresh = 0.1
print( 'LHS = ', np.abs(np.diag(matR1)) )
print( 'RHS = ', np.sum(np.abs(matR1),0) )
print( 'rat = ', np.abs(np.diag(matR1)) / np.sum(np.abs(matR1),0) )
rvcKeep = (1.0+dropThresh)*np.abs(np.diag(matR1)) > dropThresh * np.sum(np.abs(matR1),0)
print( 'rvcKeep =\n', rvcKeep )
matXK = matX[:,rvcKeep]
print( 'matX =\n', matX )
print( 'matXK =\n', matXK )
matQ, matR = linalg.qr( matXK, mode='economic' )
print( 'matQ =\n', matQ )
print( 'matR =\n', matR )
