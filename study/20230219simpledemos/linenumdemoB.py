from danutil import msg
#main_frame = inspect.currentframe()
#def msg(*arguments, **keywords):
#	print(f'[{__file__}.{inspect.stack()[1].lineno:05d}]', *arguments, **keywords)
msg('Hello from line 5 of linenumdemoB!')
def dummyfunc():
	msg('Hello from linenumdemoB.demofunc(), line 7 of linenumdemoB!')

msg('Hello from line 9 of linenumdemoB!')
