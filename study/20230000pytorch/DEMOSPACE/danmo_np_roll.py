import numpy as np

matJ = np.array([[ 11, 12, 13 ],
                 [ 21, 22, 23 ],
                 [ 31, 32, 33 ],
                 [ 41, 42, 43 ]])
print( matJ )

matJ = np.roll( matJ, 1, axis=1 )
print( matJ )
