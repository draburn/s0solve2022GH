def doubleFun( x, f ):
	return f(f(x))

def myfun( x ):
	return x**2

def returnTripleFun( f ):
	return lambda x : f(f(f(x)))

print(f'doubleFun(3.0, myfun) = {doubleFun(3.0, myfun)}')
g = returnTripleFun( myfun )
print(f'g(3.0) = {g(3.0)}')
