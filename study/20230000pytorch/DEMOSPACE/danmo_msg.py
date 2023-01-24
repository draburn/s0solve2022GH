import inspect

frame = inspect.currentframe()
def msg( *arguments, **keywords ):
	print( "[", __file__, ".", frame.f_lineno, "] ", *arguments, **keywords )

msg( "Hello world!" )
print( "__file__ = ", __file__ )
print( "frame.f_lineno = ", frame.f_lineno )
msg( "1.0+2.0+3.0 = ", 1.0+2.0+3.0, ", 1.0*2.0*3.0 = ", 1.0*2.0*3.0 )
msg( "Goodbye world!" )
