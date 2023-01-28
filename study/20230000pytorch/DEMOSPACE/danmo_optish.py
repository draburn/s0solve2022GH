def funcy( x ):
	return x**3 - 2.0
print( f'funcy(0.0) =', funcy(0.0) )
print( f'funcy(1.0) =', funcy(1.0) )
print( f'funcy(1.26) =', funcy(1.26) )
print( f'funcy(2.0) =', funcy(2.0) )

import scipy.optimize as opt
#x0, r = opt.bisect( funcy, 0.0, 2.0 )
x0 = opt.bisect( funcy, 0.0, 2.0 )
print( 'x0 = ', x0 )
#print( 'r = ', r )

exit()
# The below fails.
def funky( x ):
	return x**3
xFunky = opt.bisect( funky-2.0, 0.0, 2.0 )
print( 'xFunky = ', xFunky )
